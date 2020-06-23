#######################################
# WASH Benefits STH finished floor analysis

# distributions of EPG and Cq values
# by flooring status as box plots
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

# load data
bd <- readRDS(clean_bdata_path)
ke <- readRDS(clean_kdata_path)

bd = bd %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))
ke = ke %>% mutate(floor = factor(ifelse(floor == 1, "Improved floor", "Unimproved floor")))

#--------------------------------------
# Kato-Katz
#--------------------------------------
b_kk_plot_df = bd %>% select(floor, hwepg, ttepg) %>% mutate(country = "Bangladesh")
k_kk_plot_df = ke %>% select(floor, asca_epg) %>% mutate(country = "Kenya") %>% filter(!is.na(floor))

kk_plot_df = bind_rows(b_kk_plot_df, k_kk_plot_df)

kk_plot_df_l = melt(kk_plot_df, id.vars = c("floor", "country")) %>%
  mutate(label = case_when(
    country == "Bangladesh" & variable == "hwepg" ~ "Hookworm",
    country == "Bangladesh" & variable == "ttepg" ~ "T. trichiura",
    country == "Kenya" & variable == "asca_epg" ~ "A. lumbricoides"
  )) %>%
  filter(!is.na(value)) %>%
  mutate(label = factor(label, levels = c("T. trichiura",
                                          "Hookworm",
                                          "A. lumbricoides")))

kk_plot_b = ggplot(kk_plot_df_l %>% filter(country=="Bangladesh"), 
                 aes(x = label, y = log10(value))) +
  geom_boxplot(aes(fill = floor), width = 0.5, position=position_dodge(0.6)) +
  scale_fill_manual("", values = c("#D2D4D6", "#B06D49")) + 
  theme_bw() + xlab("") + ylab("") +
  theme(legend.position = "none") +
  scale_x_discrete(labels=make_italic) +
  scale_y_continuous(limits = c(1,5.25), labels = seq(1,5.25,0.5), breaks = seq(1,5.25,0.5)) + 
  facet_wrap(~country) + 
  theme(plot.margin = unit(c(0.2,0.2,0,0.8), "cm")) +
  coord_flip()
kk_plot_b

kk_plot_k = ggplot(kk_plot_df_l %>% filter(country=="Kenya"), 
                   aes(x = label, y = log10(value))) +
  geom_boxplot(aes(fill = floor), width = 0.5, position=position_dodge(0.6)) +
  scale_fill_manual("", values = c("#D2D4D6", "#B06D49")) + 
  theme_bw() + xlab("") + ylab("log10 eggs per gram") +
  theme(legend.position = "bottom") +
  scale_x_discrete(labels=make_italic) +
  scale_y_continuous(limits = c(1,5.25), labels = seq(1,5.25,0.5), breaks = seq(1,5.25,0.5)) + 
  facet_wrap(~country) + 
  coord_flip()
kk_plot_k

kk_plot = grid.arrange(kk_plot_b, kk_plot_k, ncol = 1, 
                      heights = c(4.8,4.3))

ggsave(filename = paste0(fig_path, "/plot_box_plot_kk.png"), plot = kk_plot,
       width = 5, height = 5)



