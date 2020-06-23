#######################################
# WASH Benefits STH finished floor analysis - Kenya

# calculate prev, geomean, CIs for KK and qPCR results
#######################################

rm(list=ls())

# Configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# Load data
ke <- readRDS(clean_kdata_path)
ke = ke %>% filter(!is.na(floor))
ke = ke %>% filter(floor != 9)

#-------------------------------------------------------------------------------
#   Convert specific data variables for analysis
#-------------------------------------------------------------------------------

# change var "floor" from continuous to factor
ke$floor <- as.factor(ke$floor)

# Log transform epg for geomean analysis
ke = ke %>%
  mutate(logalepg = log(asca_epg+1),
         logttepg = log(tric_epg+1),
         loghwepg = log(hook_epg+1))

# Prevalence outcomes list
ke_prev_outcome_kk = c("ascaris_yn","trichuris_yn","hook_yn", "sth_yn")

ke_prev_outcome_qpcr = c("pos.Al.qpcr", "pos.Tt.qpcr", "pos.Hw.qpcr",
                      "pos.Ad.qpcr", "pos.Na.qpcr", "pos.STH.qpcr","giardia_yn") # Caveat - giardia ELISA in Kenya

# Geomean outcomes list
ke_geom_outcome_kk <- c("logalepg", "logttepg", "loghwepg") 

# CT mean outcomes list                     
ke_cq_outcome <- c("al_qpcr", "tt_qpcr", "ad_qpcr", "na_qpcr") 


f1 = ke %>% filter(floor == 1)
f0 = ke %>% filter(floor == 0)

#-------------------------------------------------------------------------------
#  Calculate prevalence
#-------------------------------------------------------------------------------
ke_prev_df = ke %>%
  dplyr::select(floor, ke_prev_outcome_kk, ke_prev_outcome_qpcr) %>%
  melt(id.vars = "floor") %>%
  rename(outcome = variable) %>%
  group_by(floor, outcome) %>%
  summarise(
    N = sum(!is.na(value)),
    prev = mean(value, na.rm=T)) %>%
  mutate(diagnostic = case_when(
    outcome %in% ke_prev_outcome_kk ~ "Kato-Katz",
    outcome %in% ke_prev_outcome_qpcr ~ "qPCR"
  ))

#-------------------------------------------------------------------------------
# Get boot strap CIs
#-------------------------------------------------------------------------------
set.seed(1234)
nboot <- 1000
d = ke
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
                     pos.Al.qpcr = mean(pos.Al.qpcr, na.rm=T),
                     pos.Na.qpcr = mean(pos.Na.qpcr, na.rm=T),
                     pos.Ad.qpcr = mean(pos.Ad.qpcr, na.rm=T),
                     pos.Hw.qpcr = mean(pos.Hw.qpcr, na.rm=T),
                     pos.Tt.qpcr = mean(pos.Tt.qpcr, na.rm=T),
                     pos.STH.qpcr = mean(pos.STH.qpcr, na.rm=T),
                     giardia_yn = mean(giardia_yn, na.rm=T),
                     ascaris_yn = mean(ascaris_yn, na.rm=T),
                     trichuris_yn = mean(trichuris_yn, na.rm=T),
                     hook_yn = mean(hook_yn, na.rm=T),
                     sth_yn = mean(sth_yn, na.rm=T)) %>%
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
ke_prev_cis_df = ke_prev_df %>% 
  dplyr::select(outcome, diagnostic, floor, everything()) %>%
  left_join(prev_bs_cis, by = c("outcome", "floor")) %>%
  arrange(diagnostic, outcome, floor)

ke_prev_cis_f <- ke_prev_cis_df %>%
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

ke_geom_results_f1_list <- lapply(ke_geom_outcome_kk,
                               function(x) washb_mean(Y=f1[[x]], 
                                                      id=f1$block, print=FALSE) %>% as.data.frame())
ke_geom_results_f0_list <- lapply(ke_geom_outcome_kk,
                               function(x) washb_mean(Y=f0[[x]], 
                                                      id=f0$block, print=FALSE) %>% as.data.frame())

names(ke_geom_results_f1_list) = ke_geom_outcome_kk
names(ke_geom_results_f0_list) = ke_geom_outcome_kk

ke_geom_results_f1_df<- bind_rows(ke_geom_results_f1_list) %>% mutate(floor = 1)
ke_geom_results_f0_df <-bind_rows(ke_geom_results_f0_list) %>% mutate(floor = 0)

ke_geom_df = bind_rows(ke_geom_results_f0_df, ke_geom_results_f1_df) %>%
  mutate(outcome = rep(ke_geom_outcome_kk,2)) %>%
  # exponentiate mean, lb, ub
  mutate(Mean = exp_minus_one(Mean),
         `Lower 95%CI` = exp_minus_one(`Lower 95%CI`),
         `Upper 95%CI` = exp_minus_one(`Upper 95%CI`)) %>%
  mutate(diagnostic = "Kato-Katz")  %>%
  dplyr::select(outcome, diagnostic, floor, N, everything()) %>%
  arrange(diagnostic, outcome, floor)

ke_geom_f <- ke_geom_df %>%
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
Cq_N = ke %>% dplyr::select(floor, ke_cq_outcome) %>%
  melt(id.vars = "floor") %>%
  rename(outcome = variable) %>%
  group_by(outcome, floor) %>%
  summarise(N = sum(!is.na(value)))
 
ke_ct_f1_list <- lapply(ke_cq_outcome, 
                     function(x) summary(f1[[x]]))
ke_ct_f0_list <- lapply(ke_cq_outcome,
                     function(x) summary(f0[[x]]))

ke_ct_f1_table <- as.data.frame(do.call(rbind, ke_ct_f1_list)) %>%
  mutate(floor = 1, 
         outcome = ke_cq_outcome)
ke_ct_f0_table <- as.data.frame(do.call(rbind, ke_ct_f0_list)) %>%
  mutate(floor = 0,
         outcome = ke_cq_outcome)

ke_ct_cis_df = bind_rows(ke_ct_f1_table, ke_ct_f0_table) %>%
  mutate(diagnostic = "qPCR") %>%
  mutate(floor = as.factor(floor)) %>%
  left_join(Cq_N, by = c("outcome", "floor")) %>%
  dplyr::select(outcome, diagnostic, floor, N, everything()) %>%
  arrange(diagnostic, outcome, floor)

ke_ct_cis_f <- ke_ct_cis_df %>%
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
  
  ke_prev_cis_df, ke_prev_cis_f,
  ke_ct_cis_df, ke_ct_cis_f,
  ke_geom_df, ke_geom_f,
  
  file = paste0(results_path, "ke_prev_mean_results.RData"))




