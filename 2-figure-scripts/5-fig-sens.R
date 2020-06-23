#######################################
# WASH Benefits STH finished floor analysis

# plot results for association
# between improved floors and STH/giardia
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

results_glm = results %>% filter(fit == "GLM") 


#######################################
# plot sensitivity analyses

results_glm = results_glm %>% mutate(label = ifelse(label == "Main", "Primary analysis", label)) %>%
  filter(outcome!="positive.Hw")

plot_sens_results = function(data, diagnostic_str, analysis_str, model, heights, panel_margin, ylimits){
  plot_data = data %>% filter(diagnostic==diagnostic_str & analysis == analysis_str)
  
  scaleFUN <- function(x) sprintf("%.2f", x)
  ylabels = c(0.125, 0.25, .5, .75, 1, 1.33,2, 4)
  
  if(model == "GLM"){
    p_b = ggplot(plot_data %>% filter(country == "Bangladesh"), aes(x = outcome_f, y = IRR)) + 
      geom_point(aes(col = label), position = position_dodge(width = 0.5) ) + 
      geom_linerange(aes(col = label, ymin = `2.5%`, ymax = `97.5%`), 
                     position = position_dodge(width = 0.5)) +
      facet_wrap(~country, scales = "free_x") +
      geom_hline(yintercept = 1) + 
      scale_y_continuous(trans = 'log10', limits = ylimits, breaks = ylabels, labels = scaleFUN) +
      scale_color_discrete(name = "") +
      scale_x_discrete(labels=make_italic) +
      theme_bw() +
      theme(plot.margin = unit(c(0.2,0.2,0,panel_margin), "cm")) +
      facet_wrap(~country) + 
      ylab("") + xlab("") + coord_flip()  +
      theme(legend.position = "none")
    
    p_k = ggplot(plot_data %>% filter(country == "Kenya"), aes(x = outcome_f, y = IRR)) + 
      geom_point(aes(col = label), position = position_dodge(width = 0.5) ) + 
      geom_linerange(aes(col = label, ymin = `2.5%`, ymax = `97.5%`), 
                     position = position_dodge(width = 0.5)) +
      facet_wrap(~country, scales = "free_x") +
      geom_hline(yintercept = 1) + 
      scale_y_continuous(trans = 'log10', limits = ylimits, breaks = ylabels, labels = scaleFUN) +
      scale_color_discrete(name = "") +
      scale_x_discrete(labels=make_italic) +
      theme_bw() +
      facet_wrap(~country) + 
      ylab("Prevalence ratio (95% CI)") + xlab("") + coord_flip()  +
      theme(legend.position = "bottom")
  }
  if(model == "TMLE"){
    p_b = ggplot(plot_data %>% filter(country=="Bangladesh"), aes(x = outcome_f, y = psi)) + 
      geom_point(aes(col = label), position = position_dodge(width = 0.5) ) + 
      geom_linerange(aes(col = label, ymin = lb, ymax = ub), 
                     position = position_dodge(width = 0.5)) +
      facet_wrap(~country, scales = "free_x") +
      geom_hline(yintercept = 1) + 
      scale_y_continuous(trans = 'log10', limits = ylimits, breaks = ylabels, labels = scaleFUN) +
      scale_color_discrete(name = "") +
      scale_x_discrete(labels=make_italic) +
      theme_bw() +
      theme(plot.margin = unit(c(0.2,0.2,0,panel_margin), "cm")) +
      facet_wrap(~country) + 
      ylab("") + xlab("") + coord_flip() +
      theme(legend.position = "none")
    
    p_k = ggplot(plot_data %>% filter(country=="Kenya"), aes(x = outcome_f, y = psi)) + 
      geom_point(aes(col = label), position = position_dodge(width = 0.5) ) + 
      geom_linerange(aes(col = label, ymin = lb, ymax = ub), 
                     position = position_dodge(width = 0.5)) +
      facet_wrap(~country, scales = "free_x") +
      geom_hline(yintercept = 1) + 
      scale_y_continuous(trans = 'log10', limits = ylimits, breaks = ylabels, labels = scaleFUN) +
      scale_color_discrete(name = "") +
      scale_x_discrete(labels=make_italic) +
      theme_bw() +
      facet_wrap(~country) + 
      ylab("Prevalence ratio (95% CI)") + xlab("") + coord_flip() +
      theme(legend.position = "bottom")
  }
  
  p = grid.arrange(p_b, p_k, ncol = 1, 
                   heights = heights)
  
  return(p)
}


plot_sens_glm_qpcr_adj = plot_sens_results(data = results_glm, diagnostic_str = "qPCR", 
                     analysis_str = "Adjusted", model = "GLM",
                     heights = c(4.5,4.7), panel_margin = 0.3, 
                     ylimits = c(0.05, 4.5))
ggsave(filename = paste0(fig_path, "/plot_main_glm_qpcr_sens.png"), 
       plot = plot_sens_glm_qpcr_adj, 
       width=6, height=6)

