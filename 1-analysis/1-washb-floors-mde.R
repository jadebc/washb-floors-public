#######################################
# WASH Benefits STH finished floor analysis

# obtain minimum detectable effects
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load libraries needed for power calc only
library(lme4)
library(ICC)
library(sjstats)

#-----------------------------------------
# load data
#-----------------------------------------
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

# any STH (qpcr)
bd$positive.STH.q <- ifelse((bd$positive.Al+bd$positive.Hw+bd$positive.Tt)>0,1,0)
table(bd$positive.STH.q)

# ke
table(ke$tr)
ke$tr <- as.factor(ke$tr)
table(ke$tr)


#--------------------------------------
# power (1-beta) for a binomial outcome
# with unequal allocation
#
# adapted from Stata's power manual
# equation 1 for the X2 test
#--------------------------------------
binpower <- function(p0=0.5,p1=p0-0.05,n0=100,n1=n0,alpha=0.05) {
  za <- qnorm(1-(alpha/2))
  pbar <- (n0*p0 + n1*p1)/(n0+n1)
  sigmap <- sqrt( pbar*(1-pbar)*(1/n0 + 1/n1) )
  sigmad <- sqrt( (p0*(1-p0))/n0 + (p1*(1-p1))/n1  )
  
  term1 <- ((p1-p0) - za*sigmap) / sigmad
  term2 <- (-(p1-p0) -za*sigmap) / sigmad
  pow <- pnorm(term1) + pnorm(term2)
  return(pow)
  
}


#-----------------------------------------
# define function to obtain MDE
#-----------------------------------------
mde <- function(n0,n1,p0,alpha,beta,deff) {
  # n0    : n obs with exposure = 0
  # n1    : n obs with exposure = 1
  # p0    : P(Y|exposure=0)
  # alpha : type I error (2-tailed)
  # beta  : type II error
  # deff  : design effect: 1 + (m-1)*icc / set to 0 if no repeated measures
  
  n0adj <- round(n0/deff)
  n1adj <- round(n1/deff)
  
  p1 <- p0
  diffp <- 0.0001
  pow <- 0.1
  while(pow < (1-beta)) {
    p1 <- p1-diffp
    pow <- binpower(p0=p0,p1=p1,n0=n0adj,n1=n1adj,alpha=0.05)
  }
  
  pow_b <- pow
  mde_d <- p1-p0
  mde_r <- p1/p0
  
  res <- c(n0,n1,n0adj,n1adj,p0,p1,pow,mde_d,mde_r)
  names(res) <- c('n0','n1','n0adj','n1adj','p0','p1','power','MDEdiff','MDEratio')
  return(res)
}


#-----------------------------------------
# define function to obtain mde using 
# estimates of ICC in our dataset for a particular outcome
# and estimates of the number exposed and unexposed 
# from our data. 

# yname: name of outcome variable, as a string
# data:  data frame 
#-----------------------------------------

get_mde = function(yname, data, country, label){
  # calculate ICC from data
  icc_data = data %>% dplyr::select(!!sym(yname), clusterid)
  icc_formula = paste0(yname, "~ 1 + (1 | clusterid)") %>% as.formula()
  fit <- lmer(icc_formula,  data = icc_data)
  icc <- suppressWarnings(icc(fit))
  if(country == "Bangladesh")   deff <- 1+(7.5-1)*icc$ICC_adjusted
  if(country == "Kenya")        deff <- 1+(8.05-1)*icc$ICC_adjusted
  
  # obtain proportion exposed and unexposed
  n1_all = data %>% filter(!is.na(!!sym(yname))) %>% 
    filter(floor == 1) %>%
    summarise(n = n()) %>% 
    pull()
  
  n0_all = data %>% filter(!is.na(!!sym(yname))) %>% 
    filter(floor == 0) %>%
    summarise(n = n()) %>% 
    pull()
  
  p0_all = data %>% filter(!is.na(!!sym(yname))) %>% 
    filter(floor == 0) %>%
    summarise(mean = mean(!!sym(yname))) %>% 
    pull()
  
  n1_control = data %>% filter(!is.na(!!sym(yname))) %>% 
    filter(tr == "Control") %>%
    filter(floor == 1) %>%
    summarise(n = n()) %>% 
    pull()
  
  n0_control = data %>% filter(!is.na(!!sym(yname))) %>% 
    filter(tr == "Control") %>%
    filter(floor == 0) %>%
    summarise(n = n()) %>% 
    pull()
  
  p0_control = data %>% filter(!is.na(!!sym(yname))) %>% 
    filter(tr == "Control") %>%
    filter(floor == 0) %>%
    summarise(mean = mean(!!sym(yname))) %>% 
    pull()
  
  # MDE pooling across all arms
  mde_all = mde(
      n0 = n0_all,
      n1 = n1_all,
      p0 = p0_all,
      alpha = 0.05, 
      beta = 0.2, 
      deff = deff
    ) 
  
  # MDE in control arm only
  mde_control = mde(
    n0 = n0_control,
    n1 = n1_control,
    p0 = p0_control,
    alpha = 0.05, 
    beta = 0.2, 
    deff = deff
  ) 
  
  mde = as.data.frame(bind_rows(mde_all, mde_control)) %>%
    mutate(arms = c("All", "Control only"),
           deff = rep(deff, 2),
           country = country,
           label = label) 
  
 return(mde)
}

#-----------------------------------------
# obtain MDE
#-----------------------------------------

# bd
bd.Al.q <- get_mde("positive.Al", bd, country = "Bangladesh", label = "Ascaris qPCR")
bd.Al.kk <- get_mde("alkk", bd, country = "Bangladesh", label = "Ascaris KK")
bd.Hw.q <- get_mde("positive.Hw", bd, country = "Bangladesh", label = "Hook qPCR")
bd.Hw.kk <- get_mde("hwkk", bd, country = "Bangladesh", label = "Hook KK")
bd.Tt.q <- get_mde("positive.Tt", bd, country = "Bangladesh", label = "Trichuris qPCR")
bd.Tt.kk <- get_mde("ttkk", bd, country = "Bangladesh", label = "Trichuris KK")
bd.gia <- get_mde("posgi", bd, country = "Bangladesh", label = "Giardia qPCR")
bd.any.q <- get_mde("positive.STH.q", bd, country = "Bangladesh", label = "Any STH qPCR")

# ke
ke.Al.q <- get_mde("pos.Al.qpcr", ke, country = "Kenya", label = "Ascaris qPCR")
ke.Al.kk <- get_mde("ascaris_yn", ke, country = "Kenya", label = "Ascaris KK")
# ke.Hw.q <- get_mde("pos.Hw.qpcr",ke, country = "Kenya", label = "Hookworm qPCR") # prevalence too low
# ke.Hw.kk <- get_mde("hook_yn", ke) # prevalence too low
# ke.Tt.q <- get_mde("pos.Tt.qpcr", ke, country = "Kenya", label = "Trichuris qPCR") # prevalence too low
# ke.Tt.kk <- get_mde("trichuris_yn", ke) # prev too low
ke.gia <- get_mde("giardia_yn", ke, country = "Kenya", label = "Giardia qPCR")
ke.any.q <- get_mde("pos.STH.qpcr", ke, country = "Kenya", label = "Any STH qPCR")

#-----------------------------------------
# export as table
#-----------------------------------------

mde.bd <- bind_rows(bd.Al.q, bd.Al.kk, bd.Hw.q, bd.Hw.kk, bd.Tt.q, bd.Tt.kk, bd.gia, bd.any.q)
mde.bd[,8:9] <- format(round(mde.bd[,8:9],2), nsmall=2)
write.csv(mde.bd, file=paste0(tab_path, "/bd.mde.sth.csv"), row.names=FALSE)

mde.k <- bind_rows(ke.Al.q, ke.Al.kk, 
              ke.gia, ke.any.q)
mde.k[,8:9] <- format(round(mde.k[,8:9],2), nsmall=2)
write.csv(mde.k, file=paste0(tab_path, "/ke.mde.sth.csv"), row.names=FALSE)




