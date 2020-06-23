#######################################
# WASH Benefits STH finished floor analysis

# E-values figure
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

d = readRDS(paste0(results_path, "evalues.RDS"))

d = d %>%
  rename(outcome = yname) %>%
  add_species_names() %>% 
  add_diagnostic() %>% 
  mutate(outcome_f = factor(outcome_f, levels = c(
    "A. ceylanicum", "A. lumbricoides", "Hookworm",  
    "N. americanus", "T. trichiura", "Any STH","G. duodenalis"
  ))) %>%
  filter(outcome!="sth" & outcome!="positive.Hw" & outcome!="almh")

table = d %>% filter(fit == "GLM" & analysis=="Adjusted" & label=="Main") %>%
  dplyr::select(country, outcome_f, diagnostic, eval, eval_ci) %>%
  arrange(country, diagnostic, outcome_f)



#--------------------------------------
# save tables
#--------------------------------------

write.csv(table, file=paste0(tab_path, "/table-evalues.csv"), row.names=FALSE)


