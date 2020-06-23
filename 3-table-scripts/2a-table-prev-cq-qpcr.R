#######################################
# WASH Benefits STH finished floor analysis 

# Table of prevalence and mean Cq value - qPCR
#######################################

rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(results_path, "bd_prev_mean_results.RData"))
load(paste0(results_path, "ke_prev_mean_results.RData"))

########################################
# Format prevalence
#######################################

bd_prev_cis_f = bd_prev_cis_f %>% mutate(country="Bangladesh")
ke_prev_cis_f = ke_prev_cis_f %>% mutate(country="Kenya")

prevtable = bind_rows(bd_prev_cis_f, ke_prev_cis_f) %>%
  filter(diagnostic=="qPCR") %>%
  filter(outcome!="positive.Ad") %>%
  add_species_names() %>%
  mutate(outcome_f = factor(outcome_f, levels = c(
    "A. lumbricoides", "A. ceylanicum", "N. americanus",
    "T. trichiura", "Any STH", "G. duodenalis"
  ))) %>%
  dplyr::select(-c(diagnostic, outcome) ) %>%
  filter(!is.na(outcome_f)) %>%
  mutate(N = format(N, big.mark = ","))

########################################
# Format Cq mean
#######################################
bd_ct_cis_f = bd_ct_cis_f %>% mutate(country = "Bangladesh")
ke_ct_cis_f = ke_ct_cis_f %>% mutate(country = "Kenya")

cqtable = bind_rows(bd_ct_cis_f, ke_ct_cis_f) %>%
  filter(outcome!="CTmean.Ad" & outcome!="ad_qpcr") %>%
  add_species_names() %>%
  mutate(outcome_f = factor(outcome_f, levels = c(
    "A. lumbricoides", "A. ceylanicum", "N. americanus",
    "T. trichiura", "Any STH", "G. duodenalis"
  ))) %>%
  dplyr::select(-c(diagnostic, outcome) ) %>%
  filter(!is.na(outcome_f)) %>%
  mutate(N = format(N, big.mark = ","))


########################################
# Combine tables
#######################################
prevtable_w = pivot_wider(prevtable, names_from = floor, values_from = c(N, results))

cqtable_w = pivot_wider(cqtable, names_from = floor, values_from = c(N, results)) %>%
  rename(N1Cq = N_1,
         results_1Cq = results_1, 
         N0Cq = N_0,
         results_0Cq = results_0)

table = left_join(prevtable_w, cqtable_w, by = c("country", "outcome_f")) %>%
  dplyr::select(country, outcome_f, 
                N_1, results_1, N1Cq, results_1Cq,
                N_0, results_0, N0Cq, results_0Cq,) %>%
  arrange(country, outcome_f) 

table$N1Cq[table$outcome_f=="Any STH"] = "--"
table$N0Cq[table$outcome_f=="Any STH"] = "--"
table$results_1Cq[table$outcome_f=="Any STH"] = "--"
table$results_0Cq[table$outcome_f=="Any STH"] = "--"

table$N1Cq[table$outcome_f=="G. duodenalis"] = "--"
table$N0Cq[table$outcome_f=="G. duodenalis"] = "--"
table$results_1Cq[table$outcome_f=="G. duodenalis"] = "--"
table$results_0Cq[table$outcome_f=="G. duodenalis"] = "--"

########################################################
# Save tables
########################################################

write.csv(table, file=paste0(tab_path, "/table-qpcr-prev-cq.csv"), row.names=FALSE)
