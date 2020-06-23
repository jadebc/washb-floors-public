# Household finished flooring and soil-transmitted helminth and Giardia infections among children in rural Bangladesh and Kenya: a prospective cohort study

## Overview
Soil-transmitted helminth and *Giardia* infections pose a huge global health burden, and infections in young children can have especially negative and long-term consequences. A small but growing body of evidence suggests that household dirt floors may be an important exposure pathway. We conducted an observational analysis of data collected prospectively within the WASH Benefits randomized controlled trials to determine whether children living in rural homes with finished floors had lower levels of STH and *Giardia* infections compared to children living in homes with unfinished floors.

## Additional Information
A pre-specified analysis plan is available on [Open Science Framework](https://osf.io/dgkw5/).

To run the scripts, the data directory for the user must be changed in **`0-config.R`**. This will allow for replication of study findings using scripts in this repository. Similar directory statement changes will be needed wherever output files are saved down (e.g., raw estimates, figures). To run all scripts required to reproduce all tables and figures in the manuscript, run the bash script **`0-run-project`**.

## Directory structure
**`0-config.R` :** configuration file that sets data directories, sources base functions, and loads required libraries

**`0-base-functions` :** R script containing general functions used across the analysis

**`1-analysis` :** folder containing analysis scripts. To rerun all scripts in this subdirectory, run the bash script `0-run-analysis.sh`.    
* `1-washb-floors-mde.R` : obtain minimum detectable effects  
* `2a-washb-prev-mean-bd.R` : calculate prevalence, geometric mean, and confidence intervals for KK and qPCR data from Bangladesh  
* `2b-washb-prev-mean-ke.R` : calculate prevalence, geometric mean, and confidence intervals for KK and qPCR data from Kenya  
* `3-washb-floors-pfloors.R` : compare predicted probability (propensity scores) of an improved floor in both countries as a function of other baseline covariates   
* `4-washb-analysis.R` : fit glm and tmle to estimate association between improved floors and STH/giardia  
* `5-washb-analysis-intensity.R` : fit glm and tmle to estimate association between improved floors and STH/giardia, using quantitative measure of infection intensity  
* `6-washb-analysis-effectmod.R` : fit glm and tmle to estimate association between improved floors and STH/giardia, assessing effect modification by pre-selected covariates  
* `7-washb-analysis-pos.R` : fit glm and tmle to estimate association between improved floors and STH/giardia, excluding extreme propensity score values  
* `8-e-values.R` : estimate e-values
* `9-demographics.R` : demographic information

**`2-figure-scripts` :** folder containing figure scripts. To rerun all scripts in this subdirectory, run the bash script 0-run-figures.sh.  
* `1-figtab-prev.R` : distributions of EPG and Cq values by flooring status as table  
* `2-fig-box-plots.R` : distributions of EPG and Cq values by flooring status as box plots  
* `3-fig-box-plots.R` : distributions of EPG values by flooring status as box plots  
* `4-fig-results-effectmod-qpcr.R` : plot results for association between improved floors and STH/giardia, stratified by effect modifiers, using qPCR models (adjusted GLM models)  
* `5-fig-sens.R` : plot results from sensitivity analyses
* `6-fig-pscores.R` : plot distributions of propensity scores  

**`3-table-scripts` :** folder containing table scripts. To rerun all scripts in this subdirectory, run the bash script 0-run-figs.sh.  
* `1a-table-characteristics-bd.R` : Table 1 summary statistics for Bangladesh  
* `1b-table-characteristics-ke.R` : Table 1 summary statistics for Kenya  
* `1c-table-characteristics-full.R` : Table 1 summary statistics for both countries  
* `2a-table-prev-cq-bd.R` : generate tables for prevalence and geometric mean for Bangladesh  
* `2b-table-prev-cq-ke.R` : generate tables for prevalence and geometric mean for Kenya  
* `3a-table-characteristics-pos-bd.R` : summary statistics, sensitivity analysis trimming by extreme propensity score values for Bangladesh data  
* `3b-table-characteristics-pos-ke.R` : summary statistics, sensitivity analysis trimming by extreme propensity score values for Kenya data  
* `3c-table-characteristics-pos-full.R` : summary statistics, sensitivity analysis trimming by extreme propensity score values for both countries  
* `4-table-prim-vs-sens-adj-unadj-pr.R` : generate table comparing adjusted and unadjusted prevalence ratios between primary analysis and positivity sensitivity analysis  
* `5-table-pr-tmle.R` : generate table of unadjusted and adjusted prevalence ratios using tmle  
* `6-table-pr-glm-kk.R` : generate table of unadjusted and adjusted prevalence ratios using glm and Kato-Katz  
* `7-table-fecr-kk.R` : generate table of fecal egg count reductions using Kato-Katz  
* `8-table-fecr-qpcr-tmle.R` : generate table of fecal egg count reductions using qPCR 
* `9-table-evalues.R` : table of e-values
* `10-table-kk-intensity.R` : table of unadjusted and adjusted prevalence ratios for moderate/heavy infection intensity using Kato-Katz  

**`4-results` :** folder containing analysis results objects.

**`5-figures` :** folder containing figure files.

**`6-tables` :** folder containing table files.


Contributors: Jade Benjamin-Chung, Yoshika Crider, Stephanie Djajadi, Anmol Seth
