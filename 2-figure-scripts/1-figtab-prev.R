#######################################
# WASH Benefits STH finished floor analysis

# distributions of EPG and Cq values
# by flooring status as box plots
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
load(paste0(results_path, "bd_prev_mean_results.RData"))
load(paste0(results_path, "ke_prev_mean_results.RData"))

bd_prev_cis_df = bd_prev_cis_df %>% mutate(country="Bangladesh")
ke_prev_cis_df = ke_prev_cis_df %>% mutate(country="Kenya")

data = bind_rows(bd_prev_cis_df, ke_prev_cis_df)

data = data %>%
  filter(outcome!="alkk") %>%
  filter(outcome!="positive.Ad") %>%
  add_species_names() %>%
  mutate(outcome_f = factor(outcome_f, levels = c(
    "G. duodenalis", "Any STH", "T. trichiura",
    "Any hookworm","N. americanus", "A. ceylanicum", "A. lumbricoides"
  ))) %>%
  mutate(prev = prev*100,
         lb = lb*100,
         ub = ub*100) 


#--------------------------------------
# process analysis results
#--------------------------------------
results = readRDS(main_results_path) %>%
  rename(outcome = yname) %>%
  add_species_names() %>%
  mutate(outcome_f = factor(outcome_f, levels = c(
    "A. lumbricoides",  "A. ceylanicum", "N. americanus",
    "T. trichiura", "Any STH",  "G. duodenalis"
  ))) %>%
  add_diagnostic() 

results = results %>% filter(fit == "GLM" & diagnostic=="qPCR" & label=="Main") 
results = results %>% filter(outcome_f!="Any hookworm")
results = results %>% filter(!is.na(outcome_f))

results$result = pt.est.ci.f(mean = results$IRR,
                             lb = results$`2.5%`, 
                             ub = results$`97.5%`,
                             digits = 2,
                             scale=1) 

results = results %>% dplyr::select(country, analysis, outcome_f, N, 
                                    analysis, result) %>%
  arrange(country, outcome_f)

results_w = pivot_wider(results, names_from = analysis, values_from = c(result))

results_w$N = format(results_w$N, big.mark = ",")

#--------------------------------------
# qPCR
#--------------------------------------
qdata = data %>% filter(diagnostic == "qPCR") %>%
  ungroup() %>%
  mutate(floor = ifelse(floor == 0, "Unfinished floor", "Finished floor")) %>%
  mutate(floor = factor(floor, levels = c( "Unfinished floor", "Finished floor"))) %>%
  filter(outcome!="positive.Hw" & outcome!="pos.Hw.qpcr" &
           outcome!="pos.Ad.qpcr") 

q_plot_b= ggplot(qdata %>% filter(country=="Bangladesh"), aes(x = outcome_f, y = prev)) +
  geom_bar(aes(fill = floor), stat="identity", width = 0.5,position = position_dodge(width=0.5), col = "black", size = 0.3) +
  geom_errorbar(aes(col = floor, ymin = lb, ymax = ub),
                position = position_dodge(width=0.5), width =0.2) + 
  scale_fill_manual("", values = c("#B06D49", "#D2D4D6")) + 
  scale_color_manual("", values = c("black", "black")) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") +
  facet_wrap(~country, ncol=1, scales="free_y") +
  theme(plot.margin = unit(c(0.5,0.2,0,0.25), "cm"),
        strip.text.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12)) +
  scale_y_continuous(limits = c(0,45), breaks = seq(0,45,5), labels =seq(0,45,5))+
  scale_x_discrete(labels=make_italic) +
  coord_flip()
q_plot_b

q_plot_k= ggplot(qdata %>% filter(country=="Kenya"), aes(x = outcome_f, y = prev)) +
  geom_bar(aes(fill = floor), stat="identity", width = 0.5,position = position_dodge(width=0.5), col = "black", size = 0.3) +
  geom_errorbar(aes(col = floor, ymin = lb, ymax = ub),
                position = position_dodge(width=0.5), width =0.2) + 
  scale_fill_manual("", values = c("#B06D49", "#D2D4D6")) + 
  scale_color_manual("", values = c("black", "black")) + 
  theme_bw() + xlab("") + ylab("Prevalence (95% CI)") +
  theme(legend.position = "bottom",
        strip.text.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12)) +
  facet_wrap(~country) +
  scale_y_continuous(limits = c(-1,45), breaks = seq(0,45,5), labels =seq(0,45,5))+
  scale_x_discrete(labels=make_italic) +
  coord_flip()
q_plot_k

q_table_b = results_w %>% 
  filter(country == "Bangladesh") %>%
  select("N", "Unadjusted", "Adjusted") %>%
  mutate(N = as.character(N))

a = tableGrob(q_table_b, theme = ttheme_minimal(),
              rows = NULL, cols = c("N", "Unadjusted PR\n(95% CI)", "Adjusted PR\n(95% CI)"))
a$heights <- unit(c(0.12, rep(0.13, nrow(a) - 1), 0.09), "npc")

bangladesh_plot = grid.arrange(q_plot_b, a, ncol = 2, heights = 5)

q_table_k = results_w %>% 
  filter(country == "Kenya") %>%
  arrange(outcome_f) %>%
  select("N", "Unadjusted", "Adjusted") 

b = tableGrob(q_table_k, theme = ttheme_minimal(),
              rows = NULL, cols = c("N", "Unadjusted PR\n(95% CI)", "Adjusted PR\n(95% CI)"))

b$heights <- unit(c(0.08, rep(0.135, nrow(b) - 1), 0.22), "npc")

kenya_plot = grid.arrange(q_plot_k, b, ncol = 2, heights = c(4.8, 4))

combined_plot = grid.arrange(q_plot_b, a, q_plot_k, b, ncol = 2, heights = c(3.2,3))

ggsave(filename = paste0(fig_path, "/plot_table_prev_qpcr.png"), plot = combined_plot,
       width = 8.22, height = 9.5)

