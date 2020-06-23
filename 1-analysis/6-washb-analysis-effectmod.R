#######################################
# WASH Benefits STH finished floor analysis

# fit glm and tmle to estimate association
# between improved floors and STH/giardia
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

#----- Bangladesh -----
bd = bd %>% mutate(
  floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")),
  implatrine = as.character(implatrine)
  ) %>%
  mutate(implatrine = ifelse(implatrine == "Improved latrine", 1, 0)) %>%
  mutate(age0_5 = ifelse(agecat == "0-5 years", 1, 0))

#----- Kenya -----
ke = ke %>%
  filter(floor!=9) %>%
  mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor"))) %>%
  mutate(age0_5 = ifelse(agecat == "0-5 years", 1, 0)) %>%
  mutate(age0_5 = ifelse(is.na(agecat), NA, age0_5)) 

#-------------------------------------------------------------------------------
# Lists of outcomes and covariates
#-------------------------------------------------------------------------------
# outcome list for Bangladesh
b_ylist = list("positive.Al", "positive.Hw", "positive.Na", "positive.Ac", "positive.Tt", "positive.STH",
               "hwkk", "ttkk", "sth", "posgi")

# covariate list for Bangladesh
bd_Ws <- c("month","aged","sex","birthord","momage","momedu","momheight","Nlt18","Ncomp","hfiacat","elec","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile","tr")

# strata list for Bangladesh
b_strat_list = list("age0_5", "dw", "implatrine", "indexchild")
bd = bd %>% mutate(
  age0_5 = as.factor(age0_5),
  dw = as.factor(dw),
  implatrine = as.factor(implatrine),
  indexchild = as.factor(indexchild)
)

# kenya doesn't have birth order
k_ylist = list("ascaris_yn", "giardia_yn", "sth_yn", "pos.Al.qpcr","pos.STH.qpcr")

ke_Ws <- c("month","childage_sth","sex","momage","momedu","momheight","Nlt18","Ncomp","HHS",
           "electricity","radio", "television", "mobile", "clock", "bicycle",
           "motorcycle", "stove", "gascook","tr")

# strata list for Kenya
ke = ke %>% 
  mutate(imp_lat_bl = as.character(imp_lat_bl)) %>%
  mutate(implatrine = ifelse(imp_lat_bl == "yes", 1, 0)) %>%
  mutate(implatrine = ifelse(is.na(imp_lat_bl), NA, implatrine)) %>%
  mutate(deworm6m = as.character(deworm6m)) %>%
  mutate(dw = ifelse(deworm6m == "yes", 1, 0)) %>%
  mutate(dw = ifelse(deworm6m=="", NA, dw))

k_strat_list = list("age0_5", "dw", "implatrine", "indexchild")
ke = ke %>% mutate(
  age0_5 = as.factor(age0_5),
  dw = as.factor(dw),
  implatrine = as.factor(implatrine),
  indexchild = as.factor(indexchild)
)

#-------------------------------------------------------------------------------
# Fit GLM 
#-------------------------------------------------------------------------------
# Bangladesh
b_glm_list = list()
for(i in 1:length(b_strat_list)){
  b_glm_list[[i]] = bind_rows(lapply(b_ylist, function(x)
    fit_glm_strat(yname = x, data = bd, country = "Bangladesh", Ws = bd_Ws,
                  strat = b_strat_list[[i]]))
  )
}

b_glm_results = bind_rows(b_glm_list)

# Kenya
k_glm_list = list()
for(i in 1:length(k_strat_list)){
  print(k_strat_list[[i]])
  k_glm_list[[i]] = bind_rows(lapply(k_ylist, function(x)
    fit_glm_strat(yname = x, data = ke, country = "Kenya", Ws = ke_Ws,
                  strat = k_strat_list[[i]]))
  )
}

k_glm_results = bind_rows(k_glm_list)

#-------------------------------------------------------------------------------
# Fit TMLE
#-------------------------------------------------------------------------------
b_tmle_list = list()
for(i in 1:length(b_strat_list)){
  print(b_strat_list[i])
  b_tmle_list[[i]] = bind_rows(lapply(b_ylist, function(x) 
    fit_tmle_strat(yname = x, data = bd, country = "Bangladesh", Ws = bd_Ws, 
                  strat = b_strat_list[[i]]))
  )
}

b_tmle_results = bind_rows(b_tmle_list)


k_tmle_list = list()
for(i in 1:length(k_strat_list)){
  print(k_strat_list[i])
  k_tmle_list[[i]] = bind_rows(lapply(k_ylist, function(x)
    fit_tmle_strat(yname = x, data = ke, country = "Kenya", Ws = ke_Ws,
                   strat = k_strat_list[[i]]))
  )
}

k_tmle_results = bind_rows(k_tmle_list)


#-------------------------------------------------------------------------------
# Combine results
#-------------------------------------------------------------------------------
results = bind_rows(
  b_glm_results, k_glm_results,
  b_tmle_results, k_tmle_results
)

saveRDS(results, file = strat_results_path)

