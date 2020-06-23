#######################################
# WASH Benefits Bangladesh
# finished floor analysis

# configure data directories
# source base functions
# load libraries
#######################################
library(here)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(grid)
library(gridExtra)
library(washb)
library(SuperLearner)
library(scales)
library(sandwich)
library(lmtest)
library(zoo)
library(gam)
library(splines)
library(glmnet) 
library(arm)
library(foreach)
library(plotly) 
library(data.table)
library(tidyverse)
library(assertthat)
library(tmle)

#--------------------------------------------
# load base functions
#--------------------------------------------
source(paste0(here::here(), "/0-base-functions.R"))

#--------------------------------------------
# define clean data path
#--------------------------------------------
clean_bdata_path = paste0(here::here(),"/0-data/sth_floor_bangladesh_data.RDS")
clean_kdata_path = paste0(here::here(),"/0-data/sth_floor_kenya_data.RDS")

#--------------------------------------------
# define results paths
#--------------------------------------------
results_path = paste0(here::here(), "/4-results/")

prev_bd_path = paste0(results_path, "bd_prev_results.RDS")
gmn_bd_path = paste0(results_path, "bd_gmn_results.RDS")
prev_ke_path = paste0(results_path, "ke_prev_results.RDS")
gmn_ke_path = paste0(results_path, "ke_gmn_results.RDS")

main_results_path = paste0(results_path, "main_results.RDS")
main_results_positivity_path = paste0(results_path, "main_results_pos.RDS") 
fecr_results_path = paste0(results_path, "fecr_results.RDS")
strat_results_path = paste0(results_path, "strat_results.RDS")

bpred_path = paste0(results_path, "pred_floor_B.RDS")
kpred_path = paste0(results_path, "pred_floor_K.RDS")

#--------------------------------------------
# define figure / table paths
#--------------------------------------------
fig_path = paste0(here::here(), "/5-figures")
tab_path = paste0(here::here(), "/6-tables")


