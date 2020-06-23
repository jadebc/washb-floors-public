#######################################
# WASH Benefits STH finished floor analysis - Bangladesh

# Calculate prev, geomean, CIs for KK and qPCR results
#######################################

rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))
source(paste0(here::here(), "/0-base-functions.R"))

# load data
bd <- readRDS(clean_bdata_path)

#-------------------------------------------------------------------------------
#   Convert specific BD data variables for analysis
#-------------------------------------------------------------------------------

# change var "floor" from continuous to factor
bd$floor <- as.factor(bd$floor)

# Prevalence outcomes list (dropping alkk because of misclassification issue)
bd_prev_outcome_kk = c("ttkk","hwkk")
 
bd_prev_outcome_qpcr = c("positive.Al", "positive.Tt", "positive.Hw",
                      "positive.Ac", "positive.Ad", "positive.Na",
                      "positive.STH", "posgi")

# Geomean outcomes list
bd_geom_outcome_kk <- c("logttepg", "loghwepg") # dropping logalepg because of misclassification issues

# CT mean outcomes list                     
bd_cq_outcome <- c("CTmean.Al", "CTmean.Tt", "CTmean.Ac", "CTmean.Ad", "CTmean.Na", "ctgi") # Note no CTmean.Hw

f1 = bd %>% filter(floor == 1)
f0 = bd %>% filter(floor == 0)

#-------------------------------------------------------------------------------
#  Calculate prevalence
#-------------------------------------------------------------------------------
bd_prev_df = bd %>%
  dplyr::select(floor, bd_prev_outcome_kk, bd_prev_outcome_qpcr) %>%
  melt(id.vars = "floor") %>%
  rename(outcome = variable) %>%
  group_by(floor, outcome) %>%
  summarise(
    N = sum(!is.na(value)),
    prev = mean(value, na.rm=T)) %>%
  mutate(diagnostic = case_when(
    outcome %in% bd_prev_outcome_kk ~ "Kato-Katz",
    outcome %in% bd_prev_outcome_qpcr ~ "qPCR"
  ))

#-------------------------------------------------------------------------------
# Get boot strap CIs
#-------------------------------------------------------------------------------
set.seed(1234)
nboot <- 1000
d = bd
d = d %>% mutate(clusterid = as.character(clusterid))

# get unique list of cluster IDs
vids <- as.character(unique(d$clusterid))

# create the bootstrap samples
vsamp <- matrix(sample(vids, size = length(vids)*nboot, replace=TRUE),
                ncol = nboot,
                nrow = length(vids))

# for each bootstrap sample, estimate prevalence by flooring status
prev_boot <- foreach(booti = 1:nboot, .combine = rbind) %do% {

  di <- left_join(data.frame(clusterid = vsamp[,booti], stringsAsFactors = FALSE), 
                  d, by="clusterid") %>%
    group_by(floor) %>%
    dplyr::summarize(nobs = n(),
              positive.Al = mean(positive.Al, na.rm=T),
              positive.Na = mean(positive.Na, na.rm=T),
              positive.Ac = mean(positive.Ac, na.rm=T),
              positive.Hw = mean(positive.Hw, na.rm=T),
              positive.Ad = mean(positive.Ad, na.rm=T),
              positive.Tt = mean(positive.Tt, na.rm=T),
              positive.STH = mean(positive.STH, na.rm=T),
              posgi = mean(posgi, na.rm=T),
              ttkk = mean(ttkk, na.rm=T),
              hwkk = mean(hwkk, na.rm=T)) %>%
    mutate(brep = booti)
  di

}

prev_bootl = prev_boot %>% 
  melt(id.vars = c("floor", "nobs", "brep")) %>%
  dplyr::rename(outcome = variable, 
         prev = value)

#-------------------------------
# compute percentile 95% CIs
# merge to the sample means
#-------------------------------
prev_bs_cis = prev_bootl %>%
  ungroup() %>%
  group_by(outcome, floor) %>%
  summarise(lb = quantile(prev, probs = 0.025),
            ub = quantile(prev, probs = 0.975))

#-------------------------------
# merge with prevalence, format output
#-------------------------------
bd_prev_cis_df = bd_prev_df %>% 
  dplyr::select(outcome, diagnostic, floor, everything()) %>%
  left_join(prev_bs_cis, by = c("outcome", "floor")) %>%
  arrange(diagnostic, outcome, floor)

bd_prev_cis_f <- bd_prev_cis_df %>%
  mutate(results = 
           ptestci.format(
             x = prev, 
             lb = lb,
             ub = ub,
             decimals = 1,             
             scale = 100)
  ) %>%
  dplyr::select(outcome, diagnostic, floor, N, results)



#-------------------------------------------------------------------------------
# Column 2: Calculate geometric mean and 95%CI by flooring status (using epg values)
#-------------------------------------------------------------------------------
# geomean is the arithmetic mean of log transformed values, exponentiated back to original scale

bd_geom_results_f1_list <- lapply(bd_geom_outcome_kk,
                               function(x) washb_mean(Y=f1[[x]], 
                                                      id=f1$clusterid, print=FALSE) %>% as.data.frame())
bd_geom_results_f0_list <- lapply(bd_geom_outcome_kk,
                               function(x) washb_mean(Y=f0[[x]], 
                                                      id=f0$clusterid, print=FALSE) %>% as.data.frame())

names(bd_geom_results_f1_list) = bd_geom_outcome_kk
names(bd_geom_results_f0_list) = bd_geom_outcome_kk

bd_geom_results_f1_df<- bind_rows(bd_geom_results_f1_list) %>% mutate(floor = 1)
bd_geom_results_f0_df <-bind_rows(bd_geom_results_f0_list) %>% mutate(floor = 0)

bd_geom_df = bind_rows(bd_geom_results_f0_df, bd_geom_results_f1_df) %>%
  mutate(outcome = rep(bd_geom_outcome_kk,2)) %>%
  # exponentiate mean, lb, ub
  mutate(Mean = exp_minus_one(Mean),
         `Lower 95%CI` = exp_minus_one(`Lower 95%CI`),
         `Upper 95%CI` = exp_minus_one(`Upper 95%CI`)) %>%
  mutate(diagnostic = "Kato-Katz")  %>%
  dplyr::select(outcome, diagnostic, floor, N, everything()) %>%
  arrange(diagnostic, outcome, floor)

bd_geom_f <- bd_geom_df %>%
  mutate(results = 
           ptestci.format(
             x = Mean, 
             lb = `Lower 95%CI`,
             ub = `Upper 95%CI`,
             decimals = 2,             
             scale = 1)
  ) %>%
  dplyr::select(outcome, diagnostic, floor, N, results)

#--------------------------------------
# Column 4: Median Cq values in positive samples and range
#--------------------------------------
Cq_N = bd %>% dplyr::select(floor, bd_cq_outcome) %>%
  melt(id.vars = "floor") %>%
  rename(outcome = variable) %>%
  group_by(outcome, floor) %>%
  summarise(N = sum(!is.na(value)))

bd_ct_f1_list <- lapply(bd_cq_outcome, 
                     function(x) summary(f1[[x]]))
bd_ct_f0_list <- lapply(bd_cq_outcome,
                     function(x) summary(f0[[x]]))

bd_ct_f1_table <- as.data.frame(do.call(rbind, bd_ct_f1_list)) %>%
  mutate(floor = 1, 
         outcome = bd_cq_outcome)
bd_ct_f0_table <- as.data.frame(do.call(rbind, bd_ct_f0_list)) %>%
  mutate(floor = 0,
         outcome = bd_cq_outcome)

bd_ct_cis_df = bind_rows(bd_ct_f1_table, bd_ct_f0_table) %>%
  mutate(diagnostic = "qPCR") %>%
  mutate(floor = as.factor(floor)) %>%
  left_join(Cq_N, by = c("outcome", "floor")) %>%
  dplyr::select(outcome, diagnostic, floor, N, everything()) %>%
  arrange(diagnostic, outcome, floor)

bd_ct_cis_f <- bd_ct_cis_df %>%
  mutate(results = 
           ptestci.format(
             x = Median, 
             lb = `Min.`,
             ub = Max.,
             decimals = 1,             
             scale = 1)
  ) %>%
  dplyr::select(outcome, diagnostic, floor, N, results)

#--------------------------------
# save output
#--------------------------------
save(
  bd_prev_cis_df, bd_prev_cis_f,
  bd_ct_cis_df, bd_ct_cis_f,
  bd_geom_df, bd_geom_f,
    
  file = paste0(results_path, "bd_prev_mean_results.RData"))


