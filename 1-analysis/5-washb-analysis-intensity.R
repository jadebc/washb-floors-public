#######################################
# WASH Benefits STH finished floor analysis

# fit glm and tmle to estimate association
# between improved floors and STH/giardia

# analysis using quantitive measures of 
# infection intensity
#######################################

rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

# -------- Bangladesh -------- 
bd = bd %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

bd = bd %>% mutate(
  hwepg = ifelse(hwkk==0, NA, hwepg),
  ttepg = ifelse(ttkk==0, NA, ttepg),
  ctgi = ifelse(posgi==0, NA, ctgi)
  )

# -------- Kenya -------- 
ke = ke %>%
  filter(floor!=9) %>%
  mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

ke = ke %>% mutate(
  asca_epg = ifelse(ascaris_yn==0, NA, asca_epg)
)

#-------------------------------------------------------------------------------
# Lists of outcomes and covariates
#-------------------------------------------------------------------------------
# outcome list for Bangladesh
b_ylist = list("CTmean.Al", "CTmean.Ac", "CTmean.Na", "CTmean.Tt",
               "hwepg", "ttepg", "ctgi")

# covariate list for Bangladesh
bd_Ws <- c("month","aged","sex","birthord","momage","momedu","momheight","Nlt18","Ncomp","hfiacat","elec","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile","tr")

# kenya doesn't have birth order
k_ylist = list("asca_epg", "al_qpcr", "na_qpcr", "tt_qpcr")

ke_Ws <- c("month","childage_sth","sex","momage","momedu","momheight","Nlt18","Ncomp","HHS",
           "electricity", "radio", "television", "mobile", "clock", "bicycle",
           "motorcycle", "stove", "gascook","tr")

#-------------------------------------------------------------------------------
# Fit GLM
#-------------------------------------------------------------------------------
# Bangladesh
b_glm_results_list = lapply(b_ylist, function(x) 
  fit_glm(yname = x, data = bd, country = "Bangladesh", FECR = TRUE,
          Ws = bd_Ws, label = "Main")
)

b_glm_results = bind_rows(b_glm_results_list)

# Kenya
k_glm_results_list = lapply(k_ylist, function(x) 
  fit_glm(yname = x, data = ke, country = "Kenya",  FECR = TRUE,
          Ws = ke_Ws, label = "Main")
)

k_glm_results = bind_rows(k_glm_results_list)


#-------------------------------------------------------------------------------
# Fit TMLE
#-------------------------------------------------------------------------------
# Bangladesh
b_tmle_results_list = lapply(b_ylist, function(x) 
    fit_tmle(yname = x, data = bd, country = "Bangladesh", FECR = TRUE,
            Ws = bd_Ws, label = "Main")
  )

b_tmle_results = bind_rows(b_tmle_results_list)

# Kenya
k_tmle_results_list = lapply(k_ylist, function(x) 
  fit_tmle(yname = x, data = ke, country = "Kenya",  FECR = TRUE,
             Ws = ke_Ws, label = "Main")
  )

k_tmle_results = bind_rows(k_tmle_results_list)


#-------------------------------------------------------------------------------
# Combine results
#-------------------------------------------------------------------------------
results = bind_rows(b_glm_results, k_glm_results, b_tmle_results, k_tmle_results)

saveRDS(results, file = fecr_results_path)

