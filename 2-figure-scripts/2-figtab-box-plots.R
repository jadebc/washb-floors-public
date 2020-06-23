#######################################
# WASH Benefits STH finished floor analysis

# distributions of EPG and Cq values
# by flooring status as box plots
# and table of relative reduction in Cq values
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

bd = bd %>% mutate(floor = factor(ifelse(floor == 1, "Finished floor", "Unfinished floor")))
ke = ke %>% mutate(floor = factor(ifelse(floor == 1, "Finished floor", "Unfinished floor")))

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

results = results %>% filter(diagnostic=="qPCR" & label=="Main" & fit=="GLM") 

results$result = pt.est.ci.f(mean = results$psi,
                             lb = results$lb, 
                             ub = results$ub,
                             digits = 2,
                             scale=1) 

results = results %>% dplyr::select(country, analysis, outcome_f, N, 
                                    analysis, result)

results_w = pivot_wider(results, names_from = c(analysis), values_from = c(result))
#change -0.00 to 0.00
results_w$Adjusted[results_w$Adjusted == "-0.00 (-0.03, 0.03)"] <- "0.00 (-0.03, 0.03)"
results_w$Adjusted[results_w$Adjusted == "0.11 (-0.00, 0.22)"] <- "0.11 (0.00, 0.22)"
results_w$Unadjusted[results_w$Unadjusted == "-0.00 (-0.07, 0.06)"] <- "0.00 (-0.07, 0.06)"
results_w$N = format(results_w$N, big.mark = ",")

#--------------------------------------
# qPCR
#--------------------------------------
b_q_plot_df = bd %>% select(floor, CTmean.Al, CTmean.Na, CTmean.Ad, CTmean.Ac, CTmean.Tt, ctgi) %>%
  mutate(country = "Bangladesh")

k_q_plot_df = ke %>% select(floor, al_qpcr, na_qpcr, tt_qpcr) %>% mutate(country = "Kenya") %>% 
  filter(!is.na(floor))

q_plot_df = bind_rows(b_q_plot_df, k_q_plot_df)

q_plot_df_l = melt(q_plot_df, id.vars = c("floor", "country")) %>%
  mutate(label = case_when(
    country == "Bangladesh" & variable == "CTmean.Al" ~ "A. lumbricoides",
    country == "Bangladesh" & variable == "CTmean.Na" ~ "N. americanus",
    country == "Bangladesh" & variable == "CTmean.Ac" ~ "A. ceylanicum",
    country == "Bangladesh" & variable == "CTmean.Tt" ~ "T. trichiura",
    country == "Bangladesh" & variable == "ctgi" ~ "G. duodenalis",
    country == "Kenya" & variable == "al_qpcr" ~ "A. lumbricoides",
    country == "Kenya" & variable == "na_qpcr" ~ "N. americanus",
    country == "Kenya" & variable == "tt_qpcr" ~ "T. trichiura"
    
  )) %>%
  mutate(label = factor(label, levels = c(
    "G. duodenalis",
    "T. trichiura",
    "N. americanus",
    "A. ceylanicum",
    "A. lumbricoides"
  ))) %>%
  filter(!is.na(value)) %>%
  filter(variable!="CTmean.Ad") %>%
  mutate(floor = factor(floor, levels = c("Unfinished floor", "Finished floor")))

  
q_plot_b = ggplot(q_plot_df_l %>% filter(country=="Bangladesh"), aes(x = label, y = value)) +
  geom_boxplot(aes(fill = floor), width = 0.5, position=position_dodge(0.6)) +
  scale_fill_manual("", values = c("#B06D49", "#D2D4D6")) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none", text = element_text(size=16)) +
  facet_wrap(~country) +
  scale_x_discrete(labels=make_italic) +
  scale_y_continuous(limits = c(15,40), labels = seq(15,40,5), breaks = seq(15,40,5)) + 
  theme(plot.margin = unit(c(0.2,0,0,0.25), "cm"),
        strip.text.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        legend.text=element_text(size=12)) +
  coord_flip()
q_plot_b

q_plot_k = ggplot(q_plot_df_l %>% filter(country=="Kenya"), aes(x = label, y = value)) +
  geom_boxplot(aes(fill = floor), width = 0.5, position=position_dodge(0.6)) +
  scale_fill_manual("", values = c("#B06D49", "#D2D4D6")) + 
  theme_bw() + xlab("") + ylab("Mean Cq value") +
  theme(legend.position = "bottom", text = element_text(size=16),
        strip.text.x = element_text(size=12),
        axis.text.y = element_text(size = 12),
        axis.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  scale_x_discrete(labels=make_italic) +
  scale_y_continuous(limits = c(15,40), labels = seq(15,40,5), breaks = seq(15,40,5)) + 
  facet_wrap(~country) +
  coord_flip()
q_plot_k

#define tableGrob theme
my_theme <- ttheme_minimal(
  core=list(fg_params=list(cex = 0.75)),
  colhead=list(fg_params=list(cex = 0.8))
)

#Bangladesh
q_table_b = results_w %>% 
  filter(country == "Bangladesh") %>%
  select("N", "Unadjusted", "Adjusted") %>%
  rbind(c("", "", ""))

bangladesh_table = tableGrob(q_table_b, theme = my_theme,
              rows = NULL, cols = c("N", "Unadjusted Cq Reduction\n(95% CI)",
                                    "Adjusted Cq Reduction\n(95% CI)"))
bangladesh_table$heights <- unit(c(0.1, rep(0.147, nrow(bangladesh_table) - 2), 0.12), "npc")

#Kenya
q_table_k = results_w %>% 
  filter(country == "Kenya") %>%
  select("N", "Unadjusted", "Adjusted") %>%
  rbind(c("", "", ""))

kenya_table = tableGrob(q_table_k, theme = my_theme,
              rows = NULL, cols = c("N", "Unadjusted Cq Reduction\n(95% CI)", "Adjusted Cq Reduction\n(95% CI)"))

kenya_table$heights <- unit(c(0.12, rep(0.17, nrow(kenya_table) - 2), 0.35), "npc")

combined_plot = grid.arrange(q_plot_b, bangladesh_table, q_plot_k, kenya_table, ncol = 2, heights = c(4.2,4))

ggsave(filename = paste0(fig_path, "/plot_table_box_plot_qpcr.png"), plot = combined_plot, width = 8, height = 6.5)
