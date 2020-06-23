#######################################
# WASH Benefits STH finished floor analysis 

# Table of prevalence and mean EPG - Kato-Katz
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
  filter(diagnostic=="Kato-Katz") %>%
  filter(outcome!="alkk") %>%
  add_species_names() %>%
  mutate(outcome_f = factor(outcome_f, levels = c(
    "A. lumbricoides", "Hookworm",
    "T. trichiura", "Any STH"
  ))) %>%
  dplyr::select(-c(diagnostic, outcome) ) %>%
  filter(!is.na(outcome_f)) %>%
  mutate(N = format(N, big.mark = ","))

########################################
# Format EPG mean
#######################################
bd_geom_f = bd_geom_f %>% mutate(country = "Bangladesh")
ke_geom_f = ke_geom_f %>% mutate(country = "Kenya") %>% 
  filter(outcome!="loghwepg" & outcome!="logttepg")

epgqtable = bind_rows(bd_geom_f, ke_geom_f) %>%
  add_species_names() %>%
  mutate(outcome_f = factor(outcome_f, levels = c(
    "A. lumbricoides", "Hookworm",
    "T. trichiura", "Any STH"
  ))) %>%
  dplyr::select(-c(diagnostic, outcome) ) %>%
  filter(!is.na(outcome_f))  %>%
  mutate(N = format(N, big.mark = ","))


########################################
# Combine tables
#######################################
prevtable_w = pivot_wider(prevtable, names_from = floor, values_from = c(N, results)) 

epgqtable_w = pivot_wider(epgqtable, names_from = floor, values_from = c(N, results)) %>%
  rename(N1epg = N_1,
         results_1epg = results_1, 
         N0epg = N_0,
         results_0epg = results_0)

table = left_join(prevtable_w, epgqtable_w, by = c("country", "outcome_f")) %>%
  dplyr::select(country, outcome_f, 
                N_1, results_1, N1epg, results_1epg,
                N_0, results_0, N0epg, results_0epg,) %>%
  arrange(country, outcome_f)

table$N1epg[table$outcome_f=="Any STH"] = "--"
table$N0epg[table$outcome_f=="Any STH"] = "--"
table$results_1epg[table$outcome_f=="Any STH"] = "--"
table$results_0epg[table$outcome_f=="Any STH"] = "--"

########################################################
# Save tables
########################################################

write.csv(table, file=paste0(tab_path, "/table-kk-prev-mean.csv"), row.names=FALSE)
