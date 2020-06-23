#######################################
# WASH Benefits STH finished floor analysis

# table of unadjusted and adjusted prevalence 
# ratios using TMLE
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------
# process analysis results
#--------------------------------------
results = readRDS(main_results_path) %>%
  rename(outcome = yname) %>%
  add_species_names() %>% 
  add_diagnostic() %>% 
  mutate(outcome_f = factor(outcome_f, levels = c(
    "G. duodenalis", "Any STH", "T. trichiura",
    "Any hookworm", "Hookworm", "N. americanus", "A. ceylanicum", "A. lumbricoides"
  ))) %>%
  filter(outcome_f!="Any hookworm")



#--------------------------------------
# format tables
#--------------------------------------
results$result = pt.est.ci.f(mean = results$psi,
                             lb = results$lb, 
                             ub = results$ub,
                             digits = 2,
                             scale=1) 


results_q = results %>% filter(fit == "TMLE" & diagnostic=="qPCR" & label=="Main") 
results_kk = results %>% filter(fit == "TMLE" & diagnostic=="Kato-Katz" & label=="Main") 


results_q = results_q %>% dplyr::select(country, analysis, outcome_f, N, 
                                    analysis, result)

results_kk = results_kk %>% dplyr::select(country, analysis, outcome_f, N, 
                                        analysis, result)


results_q_w = pivot_wider(results_q, names_from = analysis, values_from = c(result))

results_kk_w = pivot_wider(results_kk, names_from = analysis, values_from = c(result))

results_q_w$N = format(results_q_w$N, big.mark = ",")
results_kk_w$N = format(results_kk_w$N, big.mark = ",")


#--------------------------------------
# save tables
#--------------------------------------

write.csv(results_q_w, file=paste0(tab_path, "/table-PRs-tmle-qpcr.csv"), row.names=FALSE)
write.csv(results_kk_w, file=paste0(tab_path, "/table-PRs-tmle-kk.csv"), row.names=FALSE)



