#######################################
# WASH Benefits STH finished floor analysis

# table of unadjusted and adjusted prevalence 
# ratios for moderate/heavy infection 
# intensity using Kato-Katz
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

results = readRDS(main_results_path) %>%
  rename(outcome = yname) %>%
  add_species_names() %>% 
  add_diagnostic() %>% 
  mutate(outcome_f = factor(outcome_f, levels = c(
    "G. duodenalis", "Any STH", "T. trichiura", "Any hookworm",
    "Hookworm", "N. americanus", "A. ceylanicum", "A. lumbricoides"
  ))) 


saveRDS(results, file = fecr_results_path)