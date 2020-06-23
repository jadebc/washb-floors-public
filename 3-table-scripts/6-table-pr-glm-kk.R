#######################################
# WASH Benefits STH finished floor analysis

# table of unadjusted and adjusted prevalence 
# ratios using GLM and Kato-Katz
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

#--------------------------------------
# format tables
#--------------------------------------
results_glm = results %>% filter(fit == "GLM" & diagnostic=="Kato-Katz" &
                                   label=="Main") 

results_glm$result = pt.est.ci.f(mean = results_glm$IRR,
                             lb = results_glm$`2.5%`, 
                             ub = results_glm$`97.5%`,
                             digits = 2,
                             scale=1) 

results_glm = results_glm %>% dplyr::select(country, analysis, outcome_f, N, 
                                    analysis, result)

results_w = pivot_wider(results_glm, names_from = analysis, values_from = c(result))

results_w$N = format(results_w$N, big.mark = ",")


#--------------------------------------
# save tables
#--------------------------------------

write.csv(results_w, file=paste0(tab_path, "/table-PRs-glm-kk.csv"), row.names=FALSE)



