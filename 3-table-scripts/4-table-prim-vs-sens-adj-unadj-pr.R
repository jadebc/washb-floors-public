#######################################
# WASH Benefits STH finished floor analysis - Bangladesh
#
# unadj versus adj PRs from primary analysis 
# compared to positivity sensitivity analysis
#######################################

rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

#--------------------------------------
# process positivity sensitivity analysis results
#--------------------------------------

results = readRDS(main_results_positivity_path) %>%
  rename(outcome = yname) %>%
  add_species_names() %>% 
  add_diagnostic() %>% 
  mutate(outcome_f = factor(outcome_f, levels = c(
    "G. duodenalis", "Any STH", "T. trichiura", "Any hookworm",
    "Hookworm", "N. americanus", "A. ceylanicum", "A. lumbricoides"
  ))) 

results_glm_pos = results %>% filter(fit == "GLM") %>%
  filter(outcome!="sth")

#--------------------------------------
# format pos results 
#--------------------------------------

### glm
results_glm_sub_pos <- results_glm_pos %>%
  mutate(result = pt.est.ci.f(mean = results_glm_pos$IRR,
                               lb = results_glm_pos$`2.5%`, 
                               ub = results_glm_pos$`97.5%`,
                               digits = 2,
                               scale=1)
  ) %>%
  filter(label == "Main"
  ) %>%
  dplyr::select(country, outcome_f, N, result, analysis, diagnostic)

results_glm_pos_sub_w = pivot_wider(results_glm_sub_pos, names_from = analysis, values_from = c(result)) 
results_glm_pos_sub_w$N = format(results_glm_pos_sub_w$N, big.mark = ",")

# split into qpcr/KK
results_glm_pos_sub_w_kk <- results_glm_pos_sub_w[results_glm_pos_sub_w$diagnostic == "Kato-Katz",] %>%
  dplyr::select(country, outcome_f, N, Unadjusted, Adjusted)
results_glm_pos_sub_w_qpcr <- results_glm_pos_sub_w[results_glm_pos_sub_w$diagnostic == "qPCR",] %>%
  dplyr::select(country, outcome_f, N, Unadjusted, Adjusted)

#--------------------------------------
# process primary results
#--------------------------------------

results = readRDS(main_results_path) %>%
  rename(outcome = yname) %>%
  add_species_names() %>% 
  add_diagnostic() %>% 
  mutate(outcome_f = factor(outcome_f, levels = c(
    "G. duodenalis", "Any STH", "T. trichiura", "Any hookworm",
    "Hookworm", "N. americanus", "A. ceylanicum", "A. lumbricoides"
  ))) 

results_glm_primary = results %>% filter(fit == "GLM") %>%
  filter(outcome!="sth")

#--------------------------------------
# format primary results 
#--------------------------------------

### glm
results_glm_primary_sub <- results_glm_primary %>%
  mutate(result = pt.est.ci.f(mean = results_glm_primary$IRR,
                              lb = results_glm_primary$`2.5%`, 
                              ub = results_glm_primary$`97.5%`,
                              digits = 2,
                              scale=1)
  ) %>%
  filter(label == "Main"
  ) %>%
  dplyr::select(country, outcome_f, N, result, analysis, diagnostic)

results_glm_primary_sub_w = pivot_wider(results_glm_primary_sub, names_from = analysis, values_from = c(result)) 
results_glm_primary_sub_w$N = format(results_glm_primary_sub_w$N, big.mark = ",")

# split into qpcr/KK
results_glm_primary_sub_w_kk <- results_glm_primary_sub_w[results_glm_primary_sub_w$diagnostic == "Kato-Katz",] %>%
  dplyr::select(country, outcome_f, N, Unadjusted, Adjusted)
results_glm_primary_sub_w_qpcr <- results_glm_primary_sub_w[results_glm_primary_sub_w$diagnostic == "qPCR",] %>%
  dplyr::select(country, outcome_f, N, Unadjusted, Adjusted)

#--------------------------------------
# format results as table
#--------------------------------------

# qpcr table comparing primary PRs with positivity analysis PRs
qpcr_table_pri_pos <- merge(results_glm_primary_sub_w_qpcr, results_glm_pos_sub_w_qpcr,
                    by=c("country", "outcome_f"), 
                    all.y=TRUE, all.x=TRUE)

colnames(qpcr_table_pri_pos) = c("country","outcome_f","N (primary)", "Unadj (primary)", "Adj (primary)",
                                 "N (pos)", "Unadj (pos)", "Adj (pos)")
  

outcome_order <- c("A. lumbricoides", "A. ceylanicum", "N. americanus", "T. trichiura", "Any STH", "G. duodenalis")
country_order <- c("Bangladesh", "Kenya")

qpcr_table_pri_pos_ordered <- qpcr_table_pri_pos %>%
  arrange(match(outcome_f, outcome_order)) %>%
  arrange(match(country, country_order)) %>%
  drop_na(outcome_f) %>%
  filter(outcome_f != "Any hookworm")

#--------------------------------------
# save table as csv
#--------------------------------------
write.csv(qpcr_table_pri_pos_ordered, file=paste0(tab_path, "/table-pr-and-pr-pos-qpcr.csv"), row.names=FALSE)