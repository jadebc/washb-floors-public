#######################################
# WASH Benefits STH finished floor analysis

# demographic information 
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#-----------------------------------------
# load data
#-----------------------------------------
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

#-----------------------------------------
# bangladesh
#-----------------------------------------
# sth cohort
min(bd$agey)
mean(bd$agey)
max(bd$agey)

# giardia cohort
b_giardia = bd %>% filter(!is.na(posgi))
min(b_giardia$agey)
mean(b_giardia$agey)
max(b_giardia$agey)

#-----------------------------------------
# kenya
#-----------------------------------------
min(ke$childage_sth/365, na.rm=TRUE)
mean(ke$childage_sth/365, na.rm=TRUE)
max(ke$childage_sth/365, na.rm=TRUE)

#-----------------------------------------
# deworming
#-----------------------------------------
mean(bd$dw, na.rm=TRUE)
prop.table(table(ke$deworm6m, useNA='ifany'))

#-----------------------------------------
# prevalence not stratifying by flooring
#-----------------------------------------
mean(bd$positive.STH, na.rm=T)
mean(bd$posgi, na.rm=T)

mean(ke$pos.STH.qpcr, na.rm=T)
mean(ke$giardia_yn, na.rm=T)

