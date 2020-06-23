##############################################
##############################################

# Documentation: check_sparsity
# Usage: check_sparsity(data, outcome, treatment, covariates)
# Description: check a data frame for positivity violations within each 
#                 treatment, outcome, covariate combination
#
# Args/Options:   
# data:           dataset to check for sparsity
# outcome:        the outcome column name, as a string
# treatment:      the treatment column name, as a string
# covariates:     covariate column name(s), as a string

# Returns: vector of covariates (as string) that pass all checks
# Output: prints rows of the data frame with conditional prevalence < 5%
# that are dropped from returned covariates list

check_sparsity = function(data, covariates, country, tolerance = 30){
  
  data = as.data.frame(data)
  
  assert_that(all(covariates %in% colnames(data)), 
              msg = "Some covariates are not in the data.")
  
  nonbinary_b = c("tr", "quarter", "aged", "birthord", "momage", "momedu", "momheight",
                "Nlt18", "Ncomp", "hfiacat")
  
  nonbinary_k = c("tr", "quarter", "childage_sth", "momage", "momedu", "momheight",
                  "Ncomp", "HHS", "Nlt18")
  
  # drop non binary variables from covariate screening
  if(country == "Bangladesh")   covariates = covariates[-which(covariates %in% nonbinary_b)]
  
  if(country == "Kenya")   covariates = covariates[-which(covariates %in% nonbinary_k)]
  
  keep_list = list()
  for(i in 1:length(covariates)){
    x= data %>% group_by(!!sym(covariates[i]), floor) %>% dplyr::summarise(n=n())
    if(min(x$n)>tolerance){
      keep_list[[i]] = covariates[i] 
    }else{
      keep_list[[i]] = NULL
    }
  }
  
  keeps = unlist(keep_list)

  if(length(keeps) < length(covariates)){
    print(paste0("Variables that have <", 
                 tolerance,
                 "observations when stratified by floor were removed from covariates list:"))
    print(covariates[!covariates %in% keeps])
    
    if(country == "Bangladesh") return(c(keeps, nonbinary_b))
    if(country == "Kenya") return(c(keeps, nonbinary_k))
  }
  
  if(length(keeps) == length(covariates)){
    if(country == "Bangladesh") return(c(covariates, nonbinary_b))
    if(country == "Kenya") return(c(covariates, nonbinary_k))
  } 
  
}

##################################################################################
# Documentation: fit_glm
# Usage: fit_glm(yname, data, country, Ws, label)
# Description: wrapper function for covariate screening and GLM 
# Args/Options: 
#   yname:    name of column with outcome variable, as a character string
#   data:     name of data frame containing outcome, exposure, and covariates
#   country:  name of country analyzed, as a character string
#   Ws:       character vector containing covariates to be screened for adjusted model 
#   label:    label to print in output data frame 
# 
# Returns: fit for unadjusted and adjusted glm
# Output: prints fit for unadjusted and adjusted glm
##################################################################################
fit_glm = function(yname, data, country, Ws, label, FECR = F){
print(yname)
  # drop rows with missing outcome
  df = data[!is.na(data[,yname]),]
  
  df$quarter = as.factor(df$quarter)
  df$tr = as.factor(df$tr)
  
  # prescreen covariates
  assert_that(all(Ws %in% colnames(df)), msg = "some covariates missing from data")
  assert_that(is.character(Ws))
  Wdf = df[,Ws]
  if(!FECR) Wset = washb_prescreen(Y = df[,yname], Ws = Wdf, family = "binomial")
  if(FECR)  Wset = washb_prescreen(Y = df[,yname], Ws = Wdf, family = "gaussian")
  
  # unadjusted
  if(!FECR){
    fit_unadj <- washb_glm(Y=df[,yname], tr=df$floor, id=df$clusterid, 
                           family = "poisson",
                           contrast=c("Unimproved floor","Improved floor"))
  }else{
    fit_unadj <- washb_glm(Y=df[,yname], tr=df$floor, id=df$clusterid, 
                           family = "gaussian", FECR="arithmetic",
                           contrast=c("Unimproved floor","Improved floor"))
  }

  # adjusted
  if(!FECR){
    fit_adj <- washb_glm(Y=df[,yname], tr=df$floor, W=df[Wset], id=df$clusterid, 
                         family = "poisson", 
                         contrast=c("Unimproved floor","Improved floor"))
  }else{
    df = df %>% dplyr::rename(treatment = tr)
    Wset[which(Wset=="tr")] = "treatment"
    
    fit_adj <- washb_glm(Y=df[,yname], tr=df$floor, W=df[Wset], id=df$clusterid, 
                         family = "gaussian", FECR="arithmetic",
                         contrast=c("Unimproved floor","Improved floor"))
  }
  
  # results 
  unadj = fit_unadj$TR %>% mutate(
    N = nrow(df),
    analysis = "Unadjusted",
    country = country,
    yname = yname) 
  
  adj = fit_adj$TR %>% mutate(
    N = nrow(df),
    analysis = "Adjusted",
    country = country,
    yname = yname) 
  
  out = bind_rows(unadj, adj) %>% select(country, yname, analysis, N, everything()) %>%
    mutate(fit = "GLM",
           label = label)
  
  if(FECR) out = out %>% dplyr::rename(var_psi = var.psi,
                                lb = ci.lb,
                                ub = ci.ub)
  
  out = out %>% mutate(strat = "None")

  out
  
  return(out)
  
}

##################################################################################
# Documentation: fit_glm_strat
# Usage: fit_glm_strat(yname, data, country, Ws, strat)
# Description: wrapper function for covariate screening and GLM 
# Args/Options: 
#   yname:    name of column with outcome variable, as a character string
#   data:     name of data frame containing outcome, exposure, and covariates
#   country:  name of country analyzed, as a character string
#   Ws:       character vector containing covariates to be screened for adjusted model 
#   strat:    name of column with indicator variable to stratify by, as a character string
# 
# Returns: fit for unadjusted and adjusted glm
# Output: prints fit for unadjusted and adjusted glm
##################################################################################
fit_glm_strat = function(yname, data, country, Ws, strat){
  
  # drop rows with missing outcome
  df = data[!is.na(data[,yname]),]

  df$quarter = as.factor(df$quarter)
  df$tr = as.factor(df$tr)
  
  # prescreen covariates
  assert_that(all(Ws %in% colnames(df)), msg = "some covariates missing from data")
  assert_that(is.character(Ws))
  
  Wdf = df[,Ws]
  Wset = washb_prescreen(Y = df[,yname], Ws = Wdf, family = "binomial")
  
  #----------------------------------------
  # get stratified estimates
  #----------------------------------------
  # filter on stratum variable
  assert_that(setequal(unique(df[,strat][!is.na(df[,strat])]), c(0,1)))
  df1 = df[df[,strat] == 1,]
  df0 = df[df[,strat] == 0,]
  
  Wset_checked1 = check_sparsity(data = df1, 
                                 covariates = Wset,
                                 tolerance = 0.05,
                                 country = country)
  
  Wset_checked0 = check_sparsity(data = df0, 
                                 covariates = Wset,
                                 tolerance = 0.05,
                                 country = country)
  
  # unadjusted
  fit_unadj1 <- washb_glm(Y=df1[,yname], tr=df1$floor, id=df1$clusterid, family = "poisson",
                         contrast=c("Unimproved floor","Improved floor"))
  
  fit_unadj0 <- washb_glm(Y=df0[,yname], tr=df0$floor, id=df0$clusterid, family = "poisson",
                          contrast=c("Unimproved floor","Improved floor"))
  
  # adjusted
  fit_adj1 <- washb_glm(Y=df1[,yname], tr=df1$floor, W=df1[Wset_checked1], id=df1$clusterid, 
                       family = "poisson", contrast=c("Unimproved floor","Improved floor"))
  
  fit_adj0 <- washb_glm(Y=df0[,yname], tr=df0$floor, W=df0[Wset_checked0], id=df0$clusterid, 
                        family = "poisson", contrast=c("Unimproved floor","Improved floor"))
  
  #----------------------------------------
  # get interaction p-values
  #----------------------------------------
  Wset_checked = check_sparsity(data = df, 
                                covariates = Wset,
                                tolerance = 0.05,
                                country = country)
  
  df = df %>% dplyr::rename(treatment = tr)
  Wset_checked = Wset_checked[-which(Wset_checked == "tr")]
  Wset_checked = c(Wset_checked, strat, "treatment")
  
  if(strat == "age0_5" & country=="Bangladesh") Wset_checked = Wset_checked[-which(Wset_checked == "aged")]
  if(strat == "age0_5" & country=="Kenya") Wset_checked = Wset_checked[-which(Wset_checked == "childage_sth")]
  
  Wdf = df %>% dplyr::select(Wset_checked)
  
  W_unadj = data.frame(V=df[,strat])
  
  fit_unadj = washb_glm(Y=df[,yname], tr=df$floor,  W=W_unadj, id=df$clusterid,
                      family = "poisson", contrast=c("Unimproved floor","Improved floor"),
                      V = "V", FECR=NULL)
  
  fit_adj = washb_glm(Y=df[,yname], tr=df$floor, W=df[Wset_checked], id=df$clusterid, 
                  family = "poisson", contrast=c("Unimproved floor","Improved floor"),
                  V = strat, FECR=NULL)
  
  # save p-values
  coefs_unadj = summary(fit_unadj$glmModel)$coef
  intrxn_pval_unadj = coefs_unadj[grep(":V1", rownames(coefs_unadj)), "Pr(>|z|)"]
  
  coefs_adj = summary(fit_adj$glmModel)$coef 
  intrxn_pval_adj = coefs_adj[grep(":V1", rownames(coefs_adj)), "Pr(>|z|)"]
  
  #----------------------------------------
  # process output
  #----------------------------------------
  # results 
  unadj1 = fit_unadj1$TR %>% mutate(
    analysis = "Unadjusted",
    country = country,
    yname = yname,
    strat = paste0(strat, "==", 1),
    N = nrow(df1),
    intrxn_pval = intrxn_pval_unadj)
  
  unadj0 = fit_unadj0$TR %>% mutate(
    analysis = "Unadjusted",
    country = country,
    yname = yname,
    strat = paste0(strat, "==", 0),
    N = nrow(df0),
    intrxn_pval = intrxn_pval_unadj)
  
  adj1 = fit_adj1$TR %>% mutate(
    analysis = "Adjusted",
    country = country,
    yname = yname,
    strat = paste0(strat, "==", 1),
    N = nrow(df1),
    intrxn_pval = intrxn_pval_adj)
  
  adj0 = fit_adj0$TR %>% mutate(
    analysis = "Adjusted",
    country = country,
    yname = yname,
    strat = paste0(strat, "==", 0),
    N = nrow(df0),
    intrxn_pval = intrxn_pval_adj) 
  
  out = bind_rows(unadj1, unadj0, adj1, adj0) %>%
    mutate(fit = "GLM") %>% 
    select(country, yname, analysis, fit, strat, N, everything()) 
 
  
  out
  
  return(out)
  
}


##################################################################################
# Documentation: fit_tmle
# Usage: fit_tmle(yname, data, country, Ws, label)
# Description: wrapper function for covariate screening and GLM 
# Args/Options: 
#   yname:    name of column with outcome variable, as a character string
#   data:     name of data frame containing outcome, exposure, and covariates
#   country:  name of country analyzed, as a character string
#   Ws:       character vector containing covariates to be screened for adjusted model 
#   label:    label to print in output data frame 
#   FECR:     Boolean for whether to estimate fecal egg count reduction (default: FALSE)
# 
# Returns: fit for unadjusted and adjusted glm
# Output: prints fit for unadjusted and adjusted glm
##################################################################################
fit_tmle = function(yname, data, country, Ws, label, FECR = FALSE){
  
  # drop rows with missing outcome
  df = data[!is.na(data[,yname]),]
  
  df$quarter = as.factor(df$quarter)
  df$tr = as.factor(df$tr)
  
  # prescreen covariates
  assert_that(all(Ws %in% colnames(df)), msg = "some covariates missing from data")
  assert_that(is.character(Ws))
  
  Ws_checked = check_sparsity(data = df, 
                              covariates = Ws,
                             country = country)
  
  Wdf = df[,Ws_checked]
  if(!FECR) Wset = washb_prescreen(Y = df[,yname], Ws = Wdf, family = "binomial")
  if(FECR)  Wset = washb_prescreen(Y = df[,yname], Ws = Wdf, family = "gaussian")
  
  dfmod = df[,c("clusterid", yname, "floor", Wset)]
  if("FALSE" %in% names(table(complete.cases(dfmod)))){
    print(paste0(as.numeric(table(complete.cases(dfmod))["FALSE"]), 
                 " rows dropped due to missing data"))
  }
  dfmod = dfmod[complete.cases(dfmod),]
  
  # unadjusted
  if(!FECR){
    fit_unadj <- washb_tmle(Y=dfmod[,yname], tr=dfmod$floor,
                            id=dfmod$clusterid,
                            family="binomial",
                            contrast=c("Unimproved floor","Improved floor"),
                            Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  }else{
    fit_unadj <- washb_tmle(Y=dfmod[,yname], tr=dfmod$floor,
                            id=dfmod$clusterid,
                            family="gaussian",
                            FECR="arithmetic",
                            contrast=c("Unimproved floor","Improved floor"),
                            Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  }

  # adjusted
  if(!FECR){
    fit_adj <- washb_tmle(Y=dfmod[,yname], tr=dfmod$floor, W=design_matrix(dfmod[,Wset]),
                          id=dfmod$clusterid,
                          family="binomial",
                          contrast=c("Unimproved floor","Improved floor"),
                          Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  }else{
    fit_adj <- washb_tmle(Y=dfmod[,yname], tr=dfmod$floor, W=design_matrix(dfmod[,Wset]),
                          id=dfmod$clusterid,
                          family="gaussian",
                          FECR="arithmetic",
                          contrast=c("Unimproved floor","Improved floor"),
                          Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  }
  
  # results 
  if(!FECR) unadj = matrix(unlist(fit_unadj$estimates$RR[c("psi","var.log.psi","CI","pvalue")])  , 1, 5) 
  if(FECR)  unadj = matrix(unlist(fit_unadj$estimates$FECR[c("psi","var.psi","CI","pvalue")])  , 1, 5) 
  
  unadj = unadj %>%
    as.data.frame() %>%
    mutate(
      analysis = "Unadjusted",
      country = country,
      yname = yname,
      N = nrow(dfmod)
    ) %>%
    dplyr::rename(psi = V1,
           var_log_psi = V2,
           lb = V3,
           ub = V4,
           pvalue = V5)
  
  if(!FECR) adj = matrix(unlist(fit_adj$estimates$RR[c("psi","var.log.psi","CI","pvalue")])  , 1, 5) 
  if(FECR)  adj = matrix(unlist(fit_adj$estimates$FECR[c("psi","var.psi","CI","pvalue")])  , 1, 5) 
  
  adj = adj %>%
    as.data.frame() %>%
    mutate(
      analysis = "Adjusted",
      country = country,
      yname = yname,
      N = nrow(dfmod)
    ) %>%
    dplyr::rename(psi = V1,
                  var_log_psi = V2,
                  lb = V3,
                  ub = V4,
                  pvalue = V5) 
  
  out = bind_rows(unadj, adj) %>% select(country, yname, analysis, N, everything()) %>%
    mutate(fit = "TMLE")
  
  out = out %>% mutate(label = label, 
                       strat = "None")
  
  if(FECR) out = out %>% dplyr::rename(var_psi = var_log_psi)

  out
  
  return(out)
  
}

##################################################################################
# Documentation: fit_tmle_strat
# Usage: fit_tmle_strat(yname, data, country, Ws, strat)
# Description: wrapper function for covariate screening and GLM 
# Args/Options: 
#   yname:    name of column with outcome variable, as a character string
#   data:     name of data frame containing outcome, exposure, and covariates
#   country:  name of country analyzed, as a character string
#   Ws:       character vector containing covariates to be screened for adjusted model 
#   strat:    name of column with indicator variable to stratify by, as a character string
# 
# Returns: fit for unadjusted and adjusted glm
# Output: prints fit for unadjusted and adjusted glm
##################################################################################
fit_tmle_strat = function(yname, data, country, Ws, strat){
  print(yname)
  
  # drop rows with missing outcome
  df = data[!is.na(data[,yname]),]
  
  df$quarter = as.factor(df$quarter)
  df$tr = as.factor(df$tr)
  
  # prescreen covariates
  assert_that(all(Ws %in% colnames(df)), msg = "some covariates missing from data")
  assert_that(is.character(Ws))

  Wdf = df[,Ws]
  Wset = washb_prescreen(Y = df[,yname], Ws = Wdf, family = "binomial")
  
  # filter on stratum variable
  assert_that(setequal(unique(df[,strat][!is.na(df[,strat])]), c(0,1)))
  df1 = df[df[,strat] == 1 & !is.na(df[,strat]),]
  df0 = df[df[,strat] == 0 & !is.na(df[,strat]),]
  
  Wset_checked1 = check_sparsity(data = df1, 
                                 covariates = Wset,
                                 country = country)
  
  Wset_checked0 = check_sparsity(data = df0, 
                                 covariates = Wset,
                                 country = country)
  
  dfmod1 = df1[,c("clusterid", yname, "floor", Wset_checked1)]
  if("FALSE" %in% names(table(complete.cases(dfmod1)))){
    print(paste0(as.numeric(table(complete.cases(dfmod1))["FALSE"]), 
                 " rows dropped due to missing data"))
  }
  dfmod1 = dfmod1[complete.cases(dfmod1),]
  
  dfmod0 = df0[,c("clusterid", yname, "floor", Wset_checked0)]
  if("FALSE" %in% names(table(complete.cases(dfmod0)))){
    print(paste0(as.numeric(table(complete.cases(dfmod0))["FALSE"]), 
                 " rows dropped due to missing data"))
  }
  dfmod0 = dfmod0[complete.cases(dfmod0),]
  
  
  # unadjusted
  fit_unadj1 <- washb_tmle(Y=dfmod1[,yname], tr=dfmod1$floor, 
                          id=dfmod1$clusterid,
                          family="binomial",
                          contrast=c("Unimproved floor","Improved floor"),
                          Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  
  fit_unadj0 <- washb_tmle(Y=dfmod0[,yname], tr=dfmod0$floor, 
                           id=dfmod0$clusterid,
                           family="binomial",
                           contrast=c("Unimproved floor","Improved floor"),
                           Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  
  # adjusted
  fit_adj1 <- washb_tmle(Y=dfmod1[,yname], tr=dfmod1$floor,
                         W=design_matrix(dfmod1[,Wset_checked1]),
                        id=dfmod1$clusterid,
                        family="binomial",
                        contrast=c("Unimproved floor","Improved floor"),
                        Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  
  fit_adj0 <- washb_tmle(Y=dfmod0[,yname], tr=dfmod0$floor,
                         W=design_matrix(dfmod0[,Wset_checked0]),
                         id=dfmod0$clusterid,
                         family="binomial",
                         contrast=c("Unimproved floor","Improved floor"),
                         Q.SL.library=c("SL.mean","SL.glm","SL.bayesglm","SL.gam","SL.glmnet"))
  
  # results 
  unadj1 = matrix(unlist(fit_unadj1$estimates$RR[c("psi","var.log.psi","CI","pvalue")])  , 1, 5)  %>%
    as.data.frame() %>%
    mutate(
      analysis = "Unadjusted",
      country = country,
      yname = yname,
      strat = paste0(strat, "==", 1),
      N = nrow(dfmod1)
    ) %>%

    dplyr::rename(psi = V1,
           var_log_psi = V2,
           lb = V3,
           ub = V4,
           pvalue = V5)
  
  unadj0 = matrix(unlist(fit_unadj0$estimates$RR[c("psi","var.log.psi","CI","pvalue")])  , 1, 5)  %>%
    as.data.frame() %>%
    mutate(
      analysis = "Unadjusted",
      country = country,
      yname = yname,
      strat = paste0(strat, "==", 0),
      N = nrow(dfmod0)
    ) %>%

    dplyr::rename(psi = V1,
           var_log_psi = V2,
           lb = V3,
           ub = V4,
           pvalue = V5)

  adj1 = matrix(unlist(fit_adj1$estimates$RR[c("psi","var.log.psi","CI","pvalue")])  , 1, 5)  %>%
    as.data.frame() %>%
    mutate(
      analysis = "Adjusted",
      country = country,
      yname = yname,
      strat = paste0(strat, "==", 1),
      N = nrow(dfmod1)
    ) %>%

    dplyr::rename(psi = V1,
           var_log_psi = V2,
           lb = V3,
           ub = V4,
           pvalue = V5) 
  
  adj0 = matrix(unlist(fit_adj0$estimates$RR[c("psi","var.log.psi","CI","pvalue")])  , 1, 5)  %>%
    as.data.frame() %>%
    mutate(
      analysis = "Adjusted",
      country = country,
      yname = yname,
      strat = paste0(strat, "==", 0),
      N = nrow(dfmod0)
    ) %>%
    dplyr::rename(psi = V1,
           var_log_psi = V2,
           lb = V3,
           ub = V4,
           pvalue = V5) 
  
  out = bind_rows(unadj1, unadj0, adj1, adj0) %>%
    mutate(fit = "TMLE") %>%
    select(country, yname, analysis, strat, fit, N, everything()) 
  
  out
  
  return(out)
  
}



make_italic = function(x){
  ifelse(x == 'A. lumbricoides', expression(italic("A. lumbricoides")),
         ifelse(x == 'G. duodenalis',  expression(italic("G. duodenalis")),
                ifelse(x == 'T. trichiura',  expression(italic("T. trichiura")),
                       ifelse(x == 'Any hookworm',  "Any hookworm", 
                              ifelse(x == 'Hookworm',  "Hookworm", 
                                   ifelse(x == "Any STH", "Any STH",
                                     ifelse(x == "A. ceylanicum", expression(italic("A. ceylanicum")),
                                            ifelse(x == "N. americanus", expression(italic("N. americanus")), "OTHER"))))))))
}

make_italic_eval = function(x){
  
  case_when(
    x == 'G. duodenalis - qPCR' ~ expression(paste(italic('G. duodenalis'), " - qPCR")),
    x == 'T. trichiura - qPCR' ~ expression(paste(italic('T. trichiura'), " - qPCR")),
    x == 'T. trichiura - Kato-Katz' ~ expression(paste(italic('T. trichiura'), " - Kato-Katz")),
    x == 'A. lumbricoides - Kato-Katz' ~ expression(paste(italic('A. lumbricoides'), " - Kato-Katz")),
    x == 'A. lumbricoides - qPCR' ~ expression(paste(italic('A. lumbricoides'), " - qPCR")),
    x == 'A. ceylanicum - qPCR' ~ expression(paste(italic('A. ceylanicum'), " - qPCR")),
    x == 'N. americanus - qPCR' ~ expression(paste(italic('N. americanus'), " - qPCR")),
    
    x == 'Hookworm - Kato-Katz' ~ expression(paste('Hookworm - Kato-Katz')),
    x == 'Hookworm - qPCR' ~ expression(paste('Hookworm - qPCR')),
    x == 'Any STH - Kato-Katz' ~ expression(paste('Any STH - Kato-Katz')),
    x == 'Any STH - qPCR' ~ expression(paste('Any STH - qPCR'))
  )
}



##############################################
##############################################
# Base function from shoo the flu github repository

mean_se=function (Y, id=NULL, Yname=NULL, print = TRUE) { # added Yname, unsure why it's dropped here and not later
  if(!is.null(id)){
    mudat <- data.frame(id = id, Y = Y)
    n.orig <- dim(mudat)[1]
    mudat <- mudat[complete.cases(mudat), ]
    n.sub <- dim(mudat)[1]
    if (n.orig > n.sub) 
      cat("\n-----------------------------------------\nDropping", 
          n.orig - n.sub, "observations\ndue to missing values in the outcome\n", 
          "Final sample size:", n.sub, "\n-----------------------------------------\n")
    fit <- glm(Y ~ 1, family = gaussian, data = mudat)
    vcovCL <- sandwichSE(fm = fit, cluster = mudat$id)
    rfit <- coeftest(fit, vcovCL)
    lb <- rfit[1, 1] - 1.96 * rfit[1, 2]
    ub <- rfit[1, 1] + 1.96 * rfit[1, 2]
    mean_ci <- matrix(c(n.sub, rfit[1, 1], sd(mudat$Y), rfit[1,2], 
                        lb, ub), nrow = 1, ncol = 6)
  }else{
    mudat <- data.frame(Y = Y)
    n.sub <- length(Y)
    fit <- glm(Y ~ 1, family = gaussian, data = mudat)
    lb <- summary(fit)$coef[1, 1] - 1.96 * summary(fit)$coef[1, 2]
    ub <- summary(fit)$coef[1, 1] + 1.96 * summary(fit)$coef[1, 2]
    mean_ci <- matrix(c(n.sub, summary(fit)$coef[1, 1], sd(mudat$Y), summary(fit)$coef[1,2], 
                        lb, ub), nrow = 1, ncol = 6)
  }
  
  colnames(mean_ci) <- c("N", "Mean", "SD", "Robust SE", "Lower 95%CI", 
                         "Upper 95%CI") 
  if (print == TRUE) 
    print(mean_ci)
  return(mean_ci)
}

sandwichSE=function (fm, cluster) 
{
  require(sandwich, quietly = TRUE)
  require(lmtest, quietly = TRUE)
  M <- length(unique(cluster))
  N <- length(cluster)
  K <- fm$rank
  dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
  uj <- apply(estfun(fm), 2, function(x) tapply(x, cluster, 
                                                sum))
  vcovCL <- dfc * sandwich(fm, meat = crossprod(uj)/N)
  return(vcovCL)
}

# get predicted means for each stratrict 
getpred=function(fit,vcov){
  
  py1=fit$coef[1]+fit$coef[2]
  py0=fit$coef[1]
  
  names(py1)=NULL
  names(py0)=NULL
  
  # se of linear combination of b0 and b1 = 
  # b0^2*s0^2 + b1^2*s1^2 + 2*b1*b0*cov(b1,b0)
  # based on the delta method
  se.1=sqrt(vcov[1,1]+vcov[2,2]+2*vcov[1,2])
  se.0=sqrt(vcov[1,1])
  
  lb1=py1-(qnorm(0.975)*se.1)
  ub1=py1+(qnorm(0.975)*se.1)
  
  lb0=py0-(qnorm(0.975)*se.0)
  ub0=py0+(qnorm(0.975)*se.0)
  
  out=as.data.frame(rbind(c(py1,lb1,ub1),c(py0,lb0,ub0)))
  colnames(out)=c("mean","lb","ub")
  out$strat=c("OUSD","WCCUSD")
  
  return(out)
}

##############################################
##############################################

# Documentation: mean_se_strat
# Usage: mean_se_strat(Y_O, Y_1, id_O, id_1, Yname, id, print = TRUE)
# Description: Wrapper to calculate mean, SE, confidence interval
#              for an outcome in each site
# Args/Options:
# Y_O:         vector containing outcome values within unfinished floor site  
# Y_1:         vector containing outcome values within finished floor site
# id_O:        vector containing hh id values within unfinished floor site
# id_1:        vector containing hh id values within finished floor site
# Yname:       name of outcome, as a string
# print:       boolean for whether to print output (Default is TRUE)
# Returns: data frame with mean, SE, confidence interval in rows for each site
# Output: none

mean_se_strat = function(Y_0, Y_1, id_0, id_1, Yname, id, print = TRUE){
  
  row_0 = mean_se(Y = Y_0,
                  Yname = Yname,
                  id = id_0)
  row_1 = mean_se(Y = Y_1,
                  Yname = Yname,
                  id = id_1)
  
  out = setDT(bind_rows(data.table(row_0), data.table(row_1)) %>% mutate(Site = c("unf", "fin"))) ## had to add data.frame or data.table to bind_rows argument
  out2 = out[, col := Yname]
  return(out2)
  
}

##############################################
##############################################

# Documentation: pt.est.ci.f
# Usage: pt.est.ci.f(mean, lb, ub, digits, scale)
# Description: Format point estimate and confidence interval nicely
# Args/Options:
# mean:        column with value of mean
# lb:          column with value of confidence interval lower bound
# ub:          column with value of confidence interval upper bound
# digits:      number of decimal points to include
# scale:       number to scale results by (e.g., 100 for percentage)

# Returns: vector of nicely formatted point estimate and confidence interval
# Output: none
pt.est.ci.f = function(mean, lb, ub, digits, scale = 1){
  format = paste0("%0.", digits, "f")
  mean_f = sprintf(format, mean*scale)
  lb_f = sprintf(format, lb*scale)
  ub_f = sprintf(format, ub*scale)
  
  out = paste0(mean_f, " (", lb_f, ", ", ub_f, ")")
  return(out)
}

##############################################
##############################################
# Needs documentation - # from 0-base-table-functions.R, wbb-sth-qpcr (slightly different from above pt.est.ci.f)
ptestci.format=function(x,lb,ub,decimals,scale){
  x=sprintf(paste0("%0.",decimals,"f"),x*scale)
  lb=sprintf(paste0("%0.",decimals,"f"),lb*scale)
  ub=sprintf(paste0("%0.",decimals,"f"),ub*scale)
  return(paste0(x," (",lb,", ",ub, ")"))
}


##############################################
##############################################
# Documentation: pt.est.f
# Usage: pt.est.ci.f(mean, digits, scale)
# Description: Format point estimate and confidence interval nicely
# Args/Options:
# mean:        column with value of mean
# lb:          column with value of confidence interval lower bound
# ub:          column with value of confidence interval upper bound
# digits:      number of decimal points to include
# scale:       number to scale results by (e.g., 100 for percentage)

# Returns: vector of nicely formatted point estimate and confidence interval
# Output: none
pt.est.f = function(mean, digits, scale = 1){
  format = paste0("%0.", digits, "f")
  mean_f = sprintf(format, mean*scale)

  out = paste0(mean_f)
  return(out)
}

##############################################
#--------------------------------
# function for exponentiating geometric mean and bounds
# exponentiate and subtract minus one
#--------------------------------
exp_minus_one = function(x){
  exp(x) - 1
}


##############################################
add_species_names = function(data){
  data %>%  mutate(outcome_f = case_when(
    outcome == "ascaris_yn" ~ "A. lumbricoides",
    outcome == "pos.Al.qpcr" ~ "A. lumbricoides",
    outcome == "positive.Al" ~ "A. lumbricoides",
    outcome == "al_qpcr" ~ "A. lumbricoides",
    outcome == "CTmean.Al" ~ "A. lumbricoides",
    outcome == "logalepg" ~ "A. lumbricoides",
    outcome == "asca_epg" ~ "A. lumbricoides",
    
    outcome == "positive.Ac" ~ "A. ceylanicum",
    outcome == "CTmean.Ac" ~ "A. ceylanicum",
    
    outcome == "positive.Na" ~ "N. americanus",
    outcome == "pos.Na.qpcr" ~ "N. americanus",
    outcome == "CTmean.Na" ~ "N. americanus",
    outcome == "na_qpcr" ~ "N. americanus",
    
    outcome == "giardia_yn" ~ "G. duodenalis",
    outcome == "ctgi" ~ "G. duodenalis",
    outcome == "posgi" ~ "G. duodenalis",
    
    outcome == "hwkk" ~ "Hookworm",
    outcome == "loghwepg" ~ "Hookworm",
    outcome == "hwepg" ~ "Hookworm",
    
    outcome == "pos.Hw.qpcr" ~ "Any hookworm",
    outcome == "positive.Hw" ~ "Any hookworm",
    
    outcome == "ttkk" ~ "T. trichiura",
    outcome == "pos.Tt.qpcr" ~ "T. trichiura",
    outcome == "positive.Tt" ~ "T. trichiura",
    outcome == "CTmean.Tt" ~ "T. trichiura",
    outcome == "tt_qpcr" ~ "T. trichiura",
    outcome == "logttepg" ~ "T. trichiura",
    outcome == "ttepg" ~ "T. trichiura",
    
    outcome == "pos.STH.qpcr" ~ "Any STH",
    outcome == "positive.STH" ~ "Any STH",
    outcome == "sth_yn" ~ "Any STH",
    
    outcome == "giardia_yn" ~ "Giardia",
    outcome == "posgi" ~ "Giardia",
    
    
  ))
}

##############################################
add_diagnostic = function(data){
  data %>% mutate(diagnostic = case_when(
    outcome == "ascaris_yn" ~ "Kato-Katz",
    outcome == "hwkk" ~ "Kato-Katz",
    outcome == "ttkk" ~ "Kato-Katz",
    outcome == "sth_yn" ~ "Kato-Katz",
    outcome == "asca_epg" ~ "Kato-Katz",
    outcome == "tt_epg" ~ "Kato-Katz",
    outcome == "hwepg" ~ "Kato-Katz",
    outcome == "ttepg" ~ "Kato-Katz",
    
    outcome == "pos.Al.qpcr" ~ "qPCR",
    outcome == "positive.Al" ~ "qPCR",
    outcome == "positive.Ac" ~ "qPCR",
    outcome == "positive.Na" ~ "qPCR",
    outcome == "giardia_yn" ~ "qPCR",
    outcome == "posgi" ~ "qPCR",
    outcome == "pos.Hw.qpcr" ~ "qPCR",
    outcome == "positive.Hw" ~ "qPCR",
    outcome == "pos.Na.qpcr" ~ "qPCR",
    outcome == "pos.Tt.qpcr" ~ "qPCR",
    outcome == "positive.Tt" ~ "qPCR",
    outcome == "pos.STH.qpcr" ~ "qPCR",
    outcome == "positive.STH" ~ "qPCR",
    
    outcome == "CTmean.Ac" ~ "qPCR",
    outcome == "CTmean.Al" ~ "qPCR",
    outcome == "CTmean.Na" ~ "qPCR",
    outcome == "CTmean.Tt" ~ "qPCR",
    outcome == "al_qpcr" ~ "qPCR",
    outcome == "na_qpcr" ~ "qPCR",
    outcome == "tt_qpcr" ~ "qPCR",
    outcome == "ctgi" ~ "qPCR",
    
  ))
}
