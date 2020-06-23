#######################################
# WASH Benefits STH finished floor analysis

# table of FECRs using KK 
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------
# process analysis results
#--------------------------------------
results = readRDS(fecr_results_path) %>%
  rename(outcome = yname) %>%
  add_species_names() %>% 
  add_diagnostic() %>% 
  mutate(outcome_f = factor(outcome_f, levels = c(
    "G. duodenalis", "T. trichiura",
    "Hookworm", "N. americanus",  "A. ceylanicum",
    "A. lumbricoides"
  ))) 

#--------------------------------------
# format tables
#--------------------------------------
results = results %>% filter(diagnostic=="qPCR" & label=="Main" & fit=="TMLE") 

results$result = pt.est.ci.f(mean = results$psi,
                             lb = results$lb, 
                             ub = results$ub,
                             digits = 2,
                             scale=1) 

results = results %>% dplyr::select(country, analysis, outcome_f, N, 
                                    analysis, result)

results_w = pivot_wider(results, names_from = analysis, values_from = c(result))

results_w$N = format(results_w$N, big.mark = ",")


#--------------------------------------
# save tables
#--------------------------------------

write.csv(results_w, file=paste0(tab_path, "/table-PRs-fecr-qpcr-tmle.csv"), row.names=FALSE)




