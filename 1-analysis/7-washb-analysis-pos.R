#######################################
# WASH Benefits STH finished floor analysis

# fit glm and tmle to estimate association
# between improved floors and STH/giardia

# excludes extreme propensity score values
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

#-------------------------------------------------------------------------------
# Create datasets for sensitivity analysis
#-------------------------------------------------------------------------------
# Bangladesh
bd_sens1 = bd %>% filter(floor_change != "Finished at baseline only")
bd_sens1 = bd_sens1 %>% filter(floor_change != "Finished at endline only")

bd_sens2 = bd %>% mutate(floor = ifelse(floor == 0 & floor_el == 1, 1, floor))

bd = bd %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

bd_sens1 = bd_sens1 %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

bd_sens2 = bd_sens2 %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

# Kenya
ke_sens1 = ke %>% filter(floor_change != "Finished at baseline only")
ke_sens1 = ke_sens1 %>% filter(floor_change != "Finished at endline only")

ke_sens2 = ke %>% mutate(floor = ifelse(floor == 0 & floor_el == 1, 1, floor))

ke_sens1 = ke_sens1 %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

ke_sens2 = ke_sens2 %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

ke = ke %>%
  filter(floor!=9) %>%
  mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))


#-------------------------------------------------------------------------------
# Lists of outcomes and covariates
#-------------------------------------------------------------------------------
# outcome list for Bangladesh
b_ylist = list("positive.Al", "positive.Ac", "positive.Na", "positive.Hw", 
               "positive.Tt", "positive.STH",
               "hwkk", "ttkk", "sth", "posgi")

# covariate list for Bangladesh
bd_Ws <- c("quarter","aged","sex","birthord","momage","momedu","momheight","Nlt18","Ncomp","hfiacat","elec","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile","tr")

# truncate data to those without potential positivity violations
bpred = readRDS(bpred_path)
bpred$dataid = as.character(bpred$dataid)

btrim_lower = min(bpred$pred, na.rm=TRUE) + 0.01
btrim_upper = max(bpred$pred, na.rm=TRUE) - 0.01

bpred_drop = bpred %>% filter(pred<btrim_lower | pred>btrim_upper) %>%
  dplyr::select(dataid) %>% pull()

nrow(bd)
bd = bd %>% filter(!dataid %in% bpred_drop)
nrow(bd)

bdata_list = list(bd, bd_sens1, bd_sens2)
names(bdata_list) = c("Main", "Sensitivity analysis 1", "Sensitivity analysis 2")


bdata_list = list(bd, bd_sens1, bd_sens2)
names(bdata_list) = c("Main", "Sensitivity analysis 1", "Sensitivity analysis 2")

# kenya doesn't have birth order
k_ylist = list("ascaris_yn", "giardia_yn", "sth_yn", 
               "pos.Al.qpcr","pos.Na.qpcr","pos.Tt.qpcr","pos.STH.qpcr")

ke_Ws <- c("quarter","childage_sth","sex","momage","momedu","momheight","Nlt18","Ncomp","HHS",
           "electricity","radio", "television", "mobile", "clock", "bicycle",
           "motorcycle", "stove", "gascook","tr")

# truncate data to those without potential positivity violations
kpred = readRDS(kpred_path)
kpred$hhid = as.character(kpred$hhid)

ktrim_lower = min(kpred$pred, na.rm=TRUE) + 0.01
ktrim_upper = max(kpred$pred, na.rm=TRUE) - 0.01

kpred_drop = kpred %>% filter(pred<ktrim_lower | pred >ktrim_upper) %>%
  dplyr::select(hhid) %>% pull()

nrow(ke)
ke = ke %>% filter(!hhid %in% kpred_drop)
nrow(ke)

kdata_list = list(ke, ke_sens1, ke_sens2)
names(kdata_list) = c("Main", "Sensitivity analysis 1", "Sensitivity analysis 2")

#-------------------------------------------------------------------------------
# Fit GLM 
#-------------------------------------------------------------------------------
# Bangladesh
b_glm_results_list_all = list()
for(i in 1:length(bdata_list)){
  print(paste0("---------------", names(bdata_list)[i], "---------------"))
  b_glm_results_list = lapply(b_ylist, function(x) 
    fit_glm(yname = x, data = bdata_list[[i]], country = "Bangladesh", 
            Ws = bd_Ws, label = names(bdata_list)[i])
  )
  b_glm_results_list_all[[i]] = bind_rows(b_glm_results_list)
  
}

b_glm_results = bind_rows(b_glm_results_list_all)

# Kenya
k_glm_results_list_all = list()
for(i in 1:length(kdata_list)){
  print(paste0("---------------", names(kdata_list)[i], "---------------"))
  k_glm_results_list = lapply(k_ylist, function(x) 
    fit_glm(yname = x, data = kdata_list[[i]], country = "Kenya", 
            Ws = ke_Ws, label = names(kdata_list)[i])
  )
  k_glm_results_list_all[[i]] = bind_rows(k_glm_results_list)
  
}

k_glm_results = bind_rows(k_glm_results_list_all)


#-------------------------------------------------------------------------------
# Fit TMLE
#-------------------------------------------------------------------------------
# Bangladesh
b_tmle_results_list_all = list()
for(i in 1:length(bdata_list)){
  print(paste0("---------------", names(bdata_list)[i], "---------------"))
  b_tmle_results_list = lapply(b_ylist, function(x) 
    fit_tmle(yname = x, data = bdata_list[[i]], country = "Bangladesh", 
            Ws = bd_Ws, label = names(bdata_list)[i])
  )
  b_tmle_results_list_all[[i]] = bind_rows(b_tmle_results_list)
  
}

b_tmle_results = bind_rows(b_tmle_results_list_all)

# Kenya
k_tmle_results_list_all = list()
for(i in 1:length(kdata_list)){
  print(paste0("---------------", names(bdata_list)[i], "---------------"))
  k_tmle_results_list = lapply(k_ylist, function(x) 
    fit_tmle(yname = x, data = kdata_list[[i]], country = "Kenya", 
             Ws = ke_Ws, label = names(kdata_list)[i])
  )
  k_tmle_results_list_all[[i]] = bind_rows(k_tmle_results_list)
  
}

k_tmle_results = bind_rows(k_tmle_results_list_all)


#-------------------------------------------------------------------------------
# Combine results
#-------------------------------------------------------------------------------
results = bind_rows(
  b_glm_results, k_glm_results,
  b_tmle_results, k_tmle_results
)

saveRDS(results, file = main_results_positivity_path)

