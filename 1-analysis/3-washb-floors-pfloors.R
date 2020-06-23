
#----------------------------------
# washb-floors-pfloors.R
#
# compare predicted probability
# of an improved floor in 
# bangladesh and kenya
# as a function of other baseline
# covariates
#----------------------------------

#----------------------------------
# preamble
#----------------------------------
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# bright color blind palette:  https://personal.sron.nl/~pault/ 
cblack <- "#000004FF"
cblue <- "#3366AA"
cteal <- "#11AA99"
cgreen <- "#66AA55"
cchartr <- "#CCCC55"
cmagent <- "#992288"
cred <- "#EE3333"
corange <- "#EEA722"
cyellow <- "#FFEE33"
cgrey <- "#777777"

#----------------------------------
# Bangladesh
#----------------------------------

bd <- readRDS(clean_bdata_path)

bd_Ws <- c("quarter","aged","sex","birthord","momage","momedu","momheight","Nlt18","Ncomp","hfiacat","elec","asset_wardrobe","asset_table","asset_chair","asset_khat","asset_chouki","asset_tv","asset_refrig","asset_bike","asset_moto","asset_sewmach","asset_mobile","tr")

# check for sparsity
Wvars_checked = check_sparsity(data = bd, covariates = bd_Ws, country = "Bangladesh")

Wvars <- bd[, colnames(bd) %in% Wvars_checked,]
Wvars <- design_matrix(Wvars)

floor <- bd$floor
id <- bd$clusterid
fitd <- data.frame(id,floor,Wvars)
fitd$dataid = bd$dataid 
fitd <- fitd[complete.cases(fitd),]

set.seed(12345)
bSLfit <- SuperLearner(Y=fitd$floor,X=fitd[names(Wvars)],id=fitd$id,
                       family='binomial',
                       SL.library=c("SL.mean","SL.glm","SL.gam","SL.glmnet",
                                    "SL.polymars","SL.xgboost"))

bpred <- bSLfit$SL.predict
bpred_df = data.frame(dataid=fitd$dataid, pred = bpred)

saveRDS(bpred_df, bpred_path)

#----------------------------------
# Kenya
#----------------------------------

kd = readRDS(clean_kdata_path)

# excluding age due to missingness
ke_Ws <- c("quarter","childage_sth","sex","momage","momedu","momheight","Nlt18","Ncomp","HHS",
           "electricity","radio", "television", "mobile", "clock", "bicycle",
           "motorcycle", "stove", "gascook","tr")

# check for sparsity
ke_Wvars_checked = check_sparsity(data = kd, covariates = ke_Ws, country = "Kenya")

kWvars <- kd[, colnames(kd) %in% ke_Wvars_checked,]

# kWvars <- kd[, colnames(kd) %in% ke_Ws,]
kWvars$tr = as.factor(kWvars$tr)
kWvars$childage_sth = as.numeric(as.character(kWvars$childage_sth))
kWvars <- kWvars[complete.cases(kWvars),]
kWvars <- design_matrix(kWvars)

kd = kd[,c("clusterid","hhid", "floor", ke_Wvars_checked)]
kd <- kd[complete.cases(kd),]
floor <- ifelse(kd$floor>=9,NA,kd$floor)
id <- kd$clusterid
hhid <- kd$hhid
kfitd <- data.frame(id,hhid,floor,kWvars)
kfitd <- kfitd[complete.cases(kfitd),]
kfitd_hh = kfitd$hhid
kfitd = kfitd %>% dplyr::select(-hhid)

set.seed(12345)
kSLfit <- SuperLearner(Y=kfitd$floor,X=kfitd[names(kWvars)],id=kfitd$id,
                       family='binomial',
                       SL.library=c("SL.mean","SL.glm","SL.gam","SL.glmnet","SL.polymars","SL.xgboost"))

kpred <- kSLfit$SL.predict
kpred_df = data.frame(hhid=kfitd_hh, pred = kpred)

saveRDS(kpred_df, kpred_path)


