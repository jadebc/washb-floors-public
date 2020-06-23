#######################################
# WASH Benefits STH finished floor analysis

# plot results for association
# between improved floors and STH/giardia

# stratify by effect modifiers 
# qPCR outcomes, adjusted GLM models
#######################################
rm(list=ls())

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))

results = readRDS(strat_results_path) %>%
  mutate(outcome = case_when(
    yname == "ascaris_yn" ~ "A. lumbricoides",
    yname == "pos.Al.qpcr" ~ "A. lumbricoides",
    yname == "positive.Al" ~ "A. lumbricoides",
    yname == "positive.Ac" ~ "A. ceylanicum",
    yname == "positive.Na" ~ "N. americanus",
    yname == "giardia_yn" ~ "G. duodenalis",
    yname == "posgi" ~ "G. duodenalis",
    yname == "hwkk" ~ "Hookworm",
    yname == "pos.Hw.qpcr" ~ "Hookworm",
    yname == "positive.Hw" ~ "Any hookworm",
    yname == "ttkk" ~ "T. trichiura",
    yname == "pos.Tt.qpcr" ~ "T. trichiura",
    yname == "positive.Tt" ~ "T. trichiura",
    yname == "sth_yn" ~ "Any STH",
    yname == "pos.STH.qpcr" ~ "Any STH",
    yname == "positive.STH" ~ "Any STH"
  )) %>%
  mutate(outcome = factor(outcome, levels = c(
    "G. duodenalis", "Any STH", "T. trichiura",
  "Hookworm","Any hookworm",  "N. americanus", "A. ceylanicum", "A. lumbricoides"
  ))) %>%
  mutate(diagnostic = case_when(
    yname == "ascaris_yn" ~ "Kato-Katz",
    yname == "pos.Al.qpcr" ~ "qPCR",
    yname == "positive.Al" ~ "qPCR",
    yname == "positive.Ac" ~ "qPCR",
    yname == "positive.Na" ~ "qPCR",
    yname == "giardia_yn" ~ "qPCR",
    yname == "posgi" ~ "qPCR",
    yname == "hwkk" ~ "Kato-Katz",
    yname == "pos.Hw.qpcr" ~ "qPCR",
    yname == "positive.Hw" ~ "qPCR",
    yname == "ttkk" ~ "Kato-Katz",
    yname == "pos.Tt.qpcr" ~ "qPCR",
    yname == "positive.Tt" ~ "qPCR",
    yname == "sth_yn" ~ "Kato-Katz",
    yname == "pos.STH.qpcr" ~ "qPCR",
    yname == "positive.STH" ~ "qPCR"
  )) %>%
  filter(!strat %in% c("agecat_t1==0", "agecat_o1==0", "agecat_c1==0")) %>%
  mutate(strat_label = case_when(
    strat == "age0_5==1" ~ "Yes",
    strat == "age0_5==0" ~ "No",
    strat == "indexchild==1" ~ "Yes",
    strat == "indexchild==0" ~ "No",
    strat == "dw==1" ~ "Yes",
    strat == "dw==0" ~ "No",
    strat == "implatrine==1" ~ "Yes",
    strat == "implatrine==0" ~ "No"
  )) %>%
  mutate(strat_cat = case_when(
    strat == "age0_5==1" ~ "Age",
    strat == "age0_5==0" ~ "Age",
    strat == "indexchild==1" ~ "Birth cohort",
    strat == "indexchild==0" ~ "Birth cohort",
    strat == "dw==1" ~ "Deworming",
    strat == "dw==0" ~ "Deworming",
    strat == "implatrine==1" ~ "Improved latrine",
    strat == "implatrine==0" ~ "Improved latrine"
  )) %>%
  mutate(strat = case_when(
    strat == "age0_5==1" ~ "Less than 5 years",
    strat == "age0_5==0" ~ "5 years or older",
    strat == "indexchild==1" ~ "Enrolled in birth cohort",
    strat == "indexchild==0" ~ "Older child",
    strat == "dw==1" ~ "Deworming in last 6 months",
    strat == "dw==0" ~ "No deworming in last 6 months",
    strat == "implatrine==1" ~ "Improved latrine",
    strat == "implatrine==0" ~ "No improved latrine"
  )) %>%
  # dropping age < 5 cutoff since it's not very 
  # different from the index child cutoff 
  filter(strat_cat != "Age")

# drop aberrant OR estimate
drops = which(results$yname=="positive.Al" & results$strat_cat=="Deworming" &
                results$country=="Bangladesh")
results[drops,"IRR"] = NA
results[drops,"2.5%"] = NA
results[drops,"97.5"] = NA

results_glm = results %>% filter(fit == "GLM" & diagnostic == "qPCR" &
                                   analysis=="Adjusted") 


#######################################
# make plots
#######################################

bplot = function(data, include_x_label, limits = NULL){
  p = ggplot(data, aes(x = outcome, y = IRR)) + 
    geom_point(aes(col = strat), position = position_dodge(width = 0.5) ) + 
    geom_linerange(aes(col = strat, ymin = `2.5%`, ymax = `97.5%`), 
                   position = position_dodge(width = 0.5)) +
    facet_grid(~strat_cat) +
    coord_flip() +
    geom_hline(yintercept = 1) + 
    theme_bw() +
    scale_color_manual(name = "", values = c("#0C5BCD", "#87B8FE")) +
    scale_shape_discrete(name = "") 
  
  
  if(!is.null(limits)){
    p = p + scale_y_continuous(limits = c(limits[1], limits[2]),
                               trans = 'log10', 
                               breaks = c(0.0625, 0.125, 0.25, .5, 1, 2, 4, 8), 
                               labels = c("0.063", "0.125", "0.25", "0.5", "1", "2", "4", "8"))
  } 
  
  if(is.null(limits)){
    p = p +
      scale_y_continuous(trans = 'log10', 
                         breaks = c(0.0625, 0.125, 0.25, .5, 1, 2, 4, 8), 
                         labels = c("0.063", "0.125", "0.25", "0.5", "1", "2", "4", "8"))
  }   
  
  p = p +
    scale_x_discrete(labels=make_italic) +
    ylab("Adjusted PR (95% CI)") +
    xlab("") +
    theme(legend.position = "bottom",
          legend.spacing.x = unit(0.01, 'cm'),
          axis.text.x = element_text(size = 7))  + 
    guides(col = guide_legend(nrow = 2))
  
  if(!include_x_label) p = p + theme(axis.text.y=element_blank(),
                                     axis.ticks.y=element_blank()) 
  

  
  return(p)
}

#-----------------------------------------
# bangladesh age plot
#-----------------------------------------
p_b_age = bplot(data = results_glm %>% filter(
  country == "Bangladesh" & strat_cat == "Birth cohort"),
  include_x_label = T) +
  theme(plot.margin = unit(c(t = 0.2, r = 0.005, b = 0.1, l = 0.2), "cm")) 

p_b_dw = bplot(data = results_glm %>% filter(
  country == "Bangladesh" & strat_cat == "Deworming"),
  include_x_label = F,
  limits = c(0.1, 6)) +
  theme(plot.margin = unit(c(t = 0.2, r = 0.005, b = 0.1, l = 0), "cm")) 

p_b_lat = bplot(data = results_glm %>% filter(
  country == "Bangladesh" & strat_cat == "Improved latrine"),
  include_x_label = F) +
  theme(plot.margin = unit(c(t = 0.2, r = 0.2, b = 0.1, l = 0), "cm")) 

p_b = grid.arrange(p_b_age, p_b_dw, p_b_lat, ncol = 3, 
                 widths = c(3,2.4, 2.4))

ggsave(filename = paste0(fig_path, "/plot_em_glm_qpcr_b.png"), plot = p_b, 
       width=9, height=4)

#-----------------------------------------
# kenya age plot
#-----------------------------------------
p_k_age = bplot(data = results_glm %>% filter(
  country == "Kenya" & strat_cat == "Birth cohort"),
  include_x_label = T) +
  theme(plot.margin = unit(c(t = 0.2, r = 0.005, b = 0.1, l = 0.2), "cm")) 

p_k_dw = bplot(data = results_glm %>% filter(
  country == "Kenya" & strat_cat == "Deworming"),
  include_x_label = F) +
  theme(plot.margin = unit(c(t = 0.2, r = 0.005, b = 0.1, l = 0), "cm")) 

p_k_lat = bplot(data = results_glm %>% filter(
  country == "Kenya" & strat_cat == "Improved latrine"),
  include_x_label = F) +
  theme(plot.margin = unit(c(t = 0.2, r = 0.2, b = 0.1, l = 0), "cm")) 

p_k = grid.arrange(p_k_age, p_k_dw, p_k_lat, ncol = 3, 
                   widths = c(3,2.4, 2.4))

ggsave(filename = paste0(fig_path, "/plot_em_glm_qpcr_k.png"), plot = p_k, 
       width=9, height=3)



