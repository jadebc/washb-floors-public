#######################################
# WASH Benefits STH finished floor analysis

# estimate e-values 
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

results = readRDS(main_results_path)
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

##################################################################################
# Documentation: get_evalue
# Usage: get_evalue(RR, lb, ub, rare, measure)
# Description: obtain e-value from measure of association and confidence interval 
# Formulas from Vanderweele et al., 2017
# Args/Options: 
#   RR:      value of measure of association 
#   lb:      value of lower bound of 95% CI of measure of association 
#   ub:      value of upper bound of 95% CI of measure of association 
#   rare:    boolean (T/F) indicating whether outcome is rare (<15%)
#   measure: character string indicating whether the measure of association is a
#            risk ratio (RR) or odds ratio (OR)
# 
# Returns: e-value of estimate and of confidence interval 
# Output: none
##################################################################################
get_evalue = function(RR, lb, ub, rare, measure){
  protective_RR = ifelse(RR < 1, T, F)
  protective_lb = ifelse(lb < 1, T, F)
  protective_ub = ifelse(ub < 1, T, F)
  
  # formula for rare outcomes (<15%)
  if((rare & measure=="OR") | measure=="RR"){
    
    # protective RR
    if(protective_RR){
      RRstar = 1/RR
      eval = RRstar + sqrt(RRstar * (RRstar - 1))
      if(protective_ub){
        ubstar = 1/ub
        eval_ci = ubstar + sqrt(ubstar *(ubstar-1))
      }else{
        eval_ci = 1
      }
    }else{
      # non protective RR 
      eval = RR + sqrt(RR * (RR-1))
      if(protective_lb){
        eval_ci = 1
      }else{
        eval_ci = lb + sqrt(lb *(lb-1))
      }
    }
  }
  
  # formula for common outcomes (>=15%)
  if(!rare & measure=="OR"){
    if(protective_RR){
      RRstar = 1/sqrt(RR)
      eval = RRstar + sqrt(RRstar * (RRstar - 1))
      if(protective_ub){
        ubstar = 1/sqrt(ub)
        eval_ci = ubstar + sqrt(ubstar *(ubstar-1))
      }else{
        eval_ci = 1
      }
    }else{
      eval = sqrt(RR) + RR + sqrt(sqrt(RR) * (sqrt(RR)-1))
      if(protective_lb){
        eval_ci = 1
      }else{
        eval_ci = sqrt(lb) + sqrt(sqrt(lb) *(sqrt(lb)-1))
      }
    }
  }
  
  out = data.frame(eval = eval, eval_ci = eval_ci)
  return(out)
  
}


#-----------------------------------------------
# identify which outcomes are rare
#-----------------------------------------------
bprev = bd %>% dplyr::select(hwkk, ttkk, sth, positive.Al, positive.Na,
                     positive.Hw, positive.Ac, positive.Ad, positive.Tt,
                     positive.STH, posgi) %>% 
  summarise(hwkk = mean(hwkk, na.rm=TRUE),
            ttkk = mean(ttkk, na.rm=TRUE),
            sth = mean(sth, na.rm=TRUE),
            posgi = mean(posgi, na.rm=TRUE),
            positive.Al = mean(positive.Al, na.rm=TRUE),
            positive.Na = mean(positive.Na, na.rm=TRUE),
            positive.Ac = mean(positive.Ac, na.rm=TRUE),
            positive.Ad = mean(positive.Ad, na.rm=TRUE),
            positive.Hw = mean(positive.Hw, na.rm=TRUE),
            positive.Tt = mean(positive.Tt, na.rm=TRUE),
            positive.STH = mean(positive.STH, na.rm=TRUE)
            ) %>% melt() %>%
  dplyr::rename(yname = variable, 
                prev = value) %>%
  mutate(rare = ifelse(prev<0.15, T, F),
         country = "Bangladesh") %>%
  dplyr::select(-prev)

# only Na is common

kprev = ke %>% dplyr::select(ascaris_yn, sth_yn, pos.Al.qpcr, 
                             pos.Na.qpcr, pos.Tt.qpcr, pos.STH.qpcr,
                             giardia_yn) %>% 
  summarise(ascaris_yn = mean(ascaris_yn, na.rm=TRUE),
            sth_yn = mean(sth_yn, na.rm=TRUE),
            pos.Al.qpcr = mean(pos.Al.qpcr, na.rm=TRUE),
            pos.Na.qpcr = mean(pos.Na.qpcr, na.rm=TRUE),
            pos.Tt.qpcr = mean(pos.Tt.qpcr, na.rm=TRUE),
            pos.STH.qpcr = mean(pos.STH.qpcr, na.rm=TRUE),
            giardia_yn = mean(giardia_yn, na.rm=TRUE)
  ) %>% melt() %>%
  dplyr::rename(yname = variable, 
                prev = value) %>%
  mutate(rare = ifelse(prev<0.15, T, F),
         country = "Kenya") %>%
  dplyr::select(-prev) 

#----------------------------------------------
# merge rare indicator onto results
#----------------------------------------------
results = results %>% mutate(
  psi = ifelse(fit=="GLM", IRR, psi),
  lb = ifelse(fit=="GLM", `2.5%`, lb),
  ub = ifelse(fit=="GLM", `97.5%`, ub)
) %>%
  dplyr::select(country, yname, analysis, fit, label, psi, lb, ub)

d = left_join(results, bprev, by = c("country", "yname"))
d = left_join(d, kprev, by = c("country", "yname"))

d = d %>% mutate(rare.x = ifelse(is.na(rare.x), rare.y, rare.x)) %>%
  dplyr::select(-rare.y) %>% dplyr::rename(rare=rare.x)

d = d %>% mutate(measure = case_when(
  fit=="GLM" ~ "RR",
  fit=="TMLE" ~ "RR"
))

#----------------------------------------------
# obtain e-values
#----------------------------------------------
d$eval = NA
d$eval_ci = NA

for(i in 1:nrow(d)){
  d$eval[i] = get_evalue(RR = d$psi[i], lb = d$lb[i], ub = d$ub[i], 
                         rare = d$rare[i], measure = d$measure[i])$eval
  d$eval_ci[i] = get_evalue(RR = d$psi[i], lb = d$lb[i], ub = d$ub[i], 
                            rare = d$rare[i], measure = d$measure[i])$eval_ci
}

saveRDS(d, paste0(results_path, "evalues.RDS"))
