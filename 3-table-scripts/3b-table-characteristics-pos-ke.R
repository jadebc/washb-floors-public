#######################################
# WASH Benefits STH finished floor analysis - Kenya

# Table 1 summary stats

# Sensitivity analysis trimming by 
# extreme propensity score values 
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))
library(dplyr)
library(data.table)

# load data
ke <- data.table((readRDS(clean_kdata_path)))

# truncate data to those without potential positivity violations
kpred = readRDS(kpred_path)
kpred$hhid = as.character(kpred$hhid)

ktrim_lower = min(kpred$pred, na.rm=TRUE) + 0.01
ktrim_upper = max(kpred$pred, na.rm=TRUE) - 0.01

kpred_drop = kpred %>% filter(pred<ktrim_lower | pred >ktrim_upper) %>%
  dplyr::select(hhid) %>% pull()

nrow(ke)
ke = ke %>% filter(!hhid %in% kpred_drop)
nrow(ke)

#---------------------------
# List and label covariates
#---------------------------

# Assign as data table
dt <- as.data.table(ke)

# Household level characteristics
#----------------------------------------
# mom's age
#----------------------------------------

row_momage = mean_se_strat(
  Y_0 = dt$momage[dt$floor == 0],
  Y_1 = dt$momage[dt$floor == 1],
  Yname = "Mother's age, years",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# mom's educ attainment (years)
#----------------------------------------
dt = dt[, edu_primary := ifelse(momedu=="Primary" | 
                                  momedu=="IncompletePrimary", 1, 0)]
dt = dt[, edu_secondary := ifelse(momedu=="AnySecondary", 1, 0)]
dt = dt[, edu_NA := ifelse(is.na(momedu), 1, 0)]
dt = dt[, edu_missing := ifelse(momedu=="", 1, 0)]

row_momedu_primary = mean_se_strat(
  Y_0 = dt$edu_primary[dt$floor == 0],
  Y_1 = dt$edu_primary[dt$floor == 1],
  Yname = "At least some primary education",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_momedu_secondary = mean_se_strat(
  Y_0 = dt$edu_secondary[dt$floor == 0],
  Y_1 = dt$edu_secondary[dt$floor == 1],
  Yname = "At least some secondary education",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_eduNA = mean_se_strat(
  Y_0 = dt$edu_NA[dt$floor == 0],
  Y_1 = dt$edu_NA[dt$floor == 1],
  Yname = "N/A education",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_edumissing = mean_se_strat(
  Y_0 = dt$edu_missing[dt$floor == 0],
  Y_1 = dt$edu_missing[dt$floor == 1],
  Yname = "Missing",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# mom's height (cm)
#----------------------------------------

row_momheight = mean_se_strat(
  Y_0 = dt$momheight[dt$floor == 0],
  Y_1 = dt$momheight[dt$floor == 1],
  Yname = "Mother's height, cm",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# N individuals in HH <=18 yrs
#----------------------------------------

row_Nind18 = mean_se_strat(
  Y_0 = dt$Nlt18[dt$floor == 0],
  Y_1 = dt$Nlt18[dt$floor == 1],
  Yname = "# individuals living in compound <=18 yrs",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# N individuals in compound
#----------------------------------------

row_Nind = mean_se_strat(
  Y_0 = dt$Ncomp[dt$floor == 0],
  Y_1 = dt$Ncomp[dt$floor == 1],
  Yname = "Total individuals living in compound",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# HFIAS category (binary secure/insecure) -- 1 = secure?
#----------------------------------------
dt = dt[, foodsecure := ifelse(HHS == 1, 1, 0)]

row_foodsecure = mean_se_strat(
  Y_0 = dt$foodsecure[dt$floor == 0], 
  Y_1 = dt$foodsecure[dt$floor == 1],
  Yname = "Food secure",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# has electricity
#----------------------------------------

row_electricity = mean_se_strat(
  Y_0 = dt$electricity[dt$floor == 0],
  Y_1 = dt$electricity[dt$floor == 1],
  Yname = "Has electricity",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# has improved wall materials
#----------------------------------------

row_walls = mean_se_strat(
  Y_0 = dt$wall[dt$floor == 0],
  Y_1 = dt$wall[dt$floor == 1],
  Yname = "Has improved wall materials",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# has improved roof material
#----------------------------------------

row_roof = mean_se_strat(
  Y_0 = dt$roof[dt$floor == 0],
  Y_1 = dt$roof[dt$floor == 1],
  Yname = "Has improved roof material",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# assets
#----------------------------------------
row_tv = mean_se_strat(
  Y_0 = dt$television[dt$floor == 0],
  Y_1 = dt$television[dt$floor == 1],
  Yname = "Owns >=1 tv",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_bike = mean_se_strat(
  Y_0 = dt$bicycle[dt$floor == 0],
  Y_1 = dt$bicycle[dt$floor == 1],
  Yname = "Owns >=1 bicycle",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_moto = mean_se_strat(
  Y_0 = dt$motorcycle[dt$floor == 0],
  Y_1 = dt$motorcycle[dt$floor == 1],
  Yname = "Owns >=1 motorcycle",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_mobile = mean_se_strat(
  Y_0 = dt$mobile[dt$floor == 0],
  Y_1 = dt$mobile[dt$floor == 1],
  Yname = "Owns >=1 mobile phone",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_stove = mean_se_strat(
  Y_0 = dt$stove[dt$floor == 0],
  Y_1 = dt$stove[dt$floor == 1],
  Yname = "Has >=1 stove",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_gascook = mean_se_strat(
  Y_0 = dt$gascook[dt$floor == 0],
  Y_1 = dt$gascook[dt$floor == 1],
  Yname = "Has gas for cooking",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# intervention arm
#----------------------------------------
dt = dt[, arm_control := ifelse(tr=="Control", 1, 0)]
dt = dt[, arm_HW := ifelse(tr=="Handwashing", 1, 0)]
dt = dt[, arm_N := ifelse(tr=="Nutrition", 1, 0)]
dt = dt[, arm_NWSH := ifelse(tr=="Nutrition + WSH", 1, 0)]
dt = dt[, arm_S := ifelse(tr=="Sanitation", 1, 0)]
dt = dt[, arm_W := ifelse(tr=="Water", 1, 0)]
dt = dt[, arm_WSH := ifelse(tr=="WSH", 1, 0)]
dt = dt[, arm_missing := ifelse(tr=="", 1, 0)]

row_arm_control = mean_se_strat(
  Y_0 = dt$arm_control[dt$floor == 0],
  Y_1 = dt$arm_control[dt$floor == 1],
  Yname = "Control",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_arm_HW = mean_se_strat(
  Y_0 = dt$arm_HW[dt$floor == 0],
  Y_1 = dt$arm_HW[dt$floor == 1],
  Yname = "Handwashing",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_arm_N = mean_se_strat(
  Y_0 = dt$arm_N[dt$floor == 0],
  Y_1 = dt$arm_N[dt$floor == 1],
  Yname = "Nutrition",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_arm_NWSH = mean_se_strat(
  Y_0 = dt$arm_NWSH[dt$floor == 0],
  Y_1 = dt$arm_NWSH[dt$floor == 1],
  Yname = "Nutrition + WSH",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_arm_S = mean_se_strat(
  Y_0 = dt$arm_S[dt$floor == 0],
  Y_1 = dt$arm_S[dt$floor == 1],
  Yname = "Sanitation",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_arm_W = mean_se_strat(
  Y_0 = dt$arm_W[dt$floor == 0],
  Y_1 = dt$arm_W[dt$floor == 1],
  Yname = "Water",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_arm_WSH = mean_se_strat(
  Y_0 = dt$arm_WSH[dt$floor == 0],
  Y_1 = dt$arm_WSH[dt$floor == 1],
  Yname = "Combined WSH",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# child mean age, years
#----------------------------------------
dt$childageyr <- dt$childage_sth/365.25

row_childageyr = mean_se_strat(
  Y_0 = dt$childageyr[dt$floor == 0],
  Y_1 = dt$childageyr[dt$floor == 1],
  Yname = "Child age, years",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

#----------------------------------------
# child age category
#----------------------------------------
# assert_that(names(table(is.na(dt$agecat)))=="FALSE") # Note 345 missing age values
dt = dt[, age_0to5 := ifelse(agecat=="0-5 years", 1, 0)]
dt = dt[, age_6to15 := ifelse(agecat =="6-15 years", 1, 0)]

row_child0to5 = mean_se_strat(
  Y_0 = dt$age_0to5[dt$floor == 0],
  Y_1 = dt$age_0to5[dt$floor == 1],
  Yname = "Child's age 0-5 years",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

# Note BD data has 6-12 years category
row_child6to15 = mean_se_strat(
  Y_0 = dt$age_6to15[dt$floor == 0],
  Y_1 = dt$age_6to15[dt$floor == 1],
  Yname = "Child's age 6-15 years",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

#----------------------------------------
# month -- move this to table 2
#----------------------------------------

#----------------------------------------
# child sex
#----------------------------------------
assert_that(names(table(is.na(dt$sex)))=="FALSE") # but note that 343 observations have missing sex
dt = dt[, sexmale := ifelse(sex=="Male", 1, 0)]
dt = dt[, sexfemale := ifelse(sex=="Female", 1, 0)]
dt = dt[, sexmissing := ifelse(sex=="missing", 1, 0)]

row_childmale = mean_se_strat(
  Y_0 = dt$sexmale[dt$floor == 0],
  Y_1 = dt$sexmale[dt$floor == 1],
  Yname = "Male, %",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

row_childfemale = mean_se_strat(
  Y_0 = dt$sexfemale[dt$floor == 0],
  Y_1 = dt$sexfemale[dt$floor == 1],
  Yname = "Female, %",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

row_childsexmissing = mean_se_strat(
  Y_0 = dt$sexmissing[dt$floor == 0],
  Y_1 = dt$sexmissing[dt$floor == 1],
  Yname = "missing, %",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

#----------------------------------------
# child type
#----------------------------------------
dt = dt[, childT1 := ifelse(child_type=="Study child (T1)", 1, 0)]
dt = dt[, childO1 := ifelse(child_type=="Older child (O1)", 1, 0)]
dt = dt[, childC2 := ifelse(child_type=="Comparison child 2 (C2)", 1, 0)]
dt = dt[, childC1 := ifelse(child_type=="Comparison child 1 (C1)", 1, 0)]
dt = dt[, child_missing := ifelse(is.na(child_type), 1, 0)]

row_childT1 = mean_se_strat(
  Y_0 = dt$childT1[dt$floor == 0],
  Y_1 = dt$childT1[dt$floor == 1],
  Yname = "Study child",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

row_childO1 = mean_se_strat(
  Y_0 = dt$childO1[dt$floor == 0],
  Y_1 = dt$childO1[dt$floor == 1],
  Yname = "Older child",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

row_childC2 = mean_se_strat(
  Y_0 = dt$childC2[dt$floor == 0],
  Y_1 = dt$childC2[dt$floor == 1],
  Yname = "Comparison child 2 (C2)",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

row_childC1 = mean_se_strat(
  Y_0 = dt$childC1[dt$floor == 0],
  Y_1 = dt$childC1[dt$floor == 1],
  Yname = "Comparison child 1 (C1)",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

row_childmissing = mean_se_strat(
  Y_0 = dt$child_missing[dt$floor == 0],
  Y_1 = dt$child_missing[dt$floor == 1],
  Yname = "Missing child type data",
  id_0 = dt$childid[dt$floor == 0],
  id_1 = dt$childid[dt$floor == 1]
)

########################################################
# Make tables
########################################################

# Table of household characteristics
table = bind_rows(
  row_momage, row_momedu_primary, row_momedu_secondary, 
  row_momheight,
  row_Nind18, row_Nind,
  row_foodsecure,
  row_electricity, row_walls, row_roof, 
  row_tv, row_bike, row_moto, row_mobile, 
  row_childageyr, row_childmale
)

noscale <- c("Mother's age, years", "Mother's height, cm", 
             "# individuals living in compound <=18 yrs", 
             "Total individuals living in compound",
             "Child age, years"
)

table <- table %>%
  mutate(results = if_else(col %in% noscale, 
      pt.est.f(
           mean = Mean, 
           digits = 1,             
           scale = 1
          ),
      pt.est.f(
        mean = Mean,
        digits = 1,
        scale = 100
      )))

Outcome_list <- c("Mother's age, years", "Mother's height, cm", 
                  "No education", "At least some primary education", "At least some secondary education",
                  "# individuals living in compound <=18 yrs", "Total individuals living in compound",
                  "Food secure",
                  "Has electricity", "Has improved wall materials", "Has improved roof material", 
                  "Owns >=1 tv", "Owns >=1 bicycle", "Owns >=1 motorcycle", "Owns >=1 mobile phone",
                  "Child age, years", "Male, %"
)


# To fix error: Each row of output must be identified by a unique combination of keys.
table <- table %>%
  dplyr::mutate(obs = row_number())

table_wide = table %>% select(col, results, Site) %>%
  spread(Site, results)

table_wide_N = table %>% dplyr::select(col, N, Site) %>%
  spread(Site, N) 

assert_that(all(table_wide$Outcome == table_wide_N$Outcome))

table_wide_all=data.frame(cbind(
  as.character(table_wide$col),
  table_wide_N$fin, 
  table_wide$fin, 
  table_wide_N$unf, 
  table_wide$unf
))

colnames(table_wide_all) = c(
  "Variable",
  "N, finished", "Result, finished",
  "N, unfinished", "Results, unfinished"
)

# Reorder the rows
table_wide_all_ordered <- table_wide_all %>%
  arrange(match(Variable, Outcome_list))

maternal_header = c("Maternal", rep("", 4))
edu_header = c("Mother's education", rep("", 4))
compound_header = c("Compound", rep("", 4))
household_header = c("Household", rep("", 4))
arm_header = c("Intervention assignment", rep("", 4))
child_header = c("Child", rep("", 4))

names(maternal_header) = colnames(table_wide_all_ordered)
names(edu_header) = colnames(table_wide_all_ordered)
names(compound_header) = colnames(table_wide_all_ordered)
names(household_header) = colnames(table_wide_all_ordered)
names(arm_header) = colnames(table_wide_all_ordered)
names(child_header) = colnames(table_wide_all_ordered)

table_wide_out = bind_rows(
  maternal_header, table_wide_all_ordered[1:4,],
  compound_header, table_wide_all_ordered[5:6,],
  household_header, table_wide_all_ordered[7:16,],
  child_header, table_wide_all_ordered[15:16,]
)

table_ke_full <- table_wide_out

########################################################
# Save tables
########################################################

write.csv(table_ke_full, file=paste0(tab_path, "/table-characteristics-pos-ke.txt"), row.names=FALSE)
write.csv(table_ke_full, file=paste0(tab_path, "/table-characteristics-pos-ke.csv"), row.names=FALSE)

