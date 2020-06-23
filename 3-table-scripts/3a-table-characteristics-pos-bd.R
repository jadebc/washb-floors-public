#######################################
# WASH Benefits STH finished floor analysis - Bangladesh

# Table 1 summary stats

# Sensitivity analysis trimming by 
# extreme propensity score values 
#######################################

# configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))
library(dplyr)
library(data.table)

# load data
bd <- data.table((readRDS(clean_bdata_path)))

# truncate data to those without potential positivity violations
bpred = readRDS(bpred_path)
bpred$dataid = as.character(bpred$dataid)

btrim_lower = min(bpred$pred, na.rm=TRUE) + 0.05
btrim_upper = max(bpred$pred, na.rm=TRUE) - 0.05

bpred_drop = bpred %>% filter(pred<btrim_lower | pred>btrim_upper) %>%
  dplyr::select(dataid) %>% pull()

nrow(bd)
bd = bd %>% filter(!dataid %in% bpred_drop)
nrow(bd)

#---------------------------
# List and label covariates
#---------------------------

# Assign as data table
dt <- as.data.table(bd)

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
dt = dt[, edu_none := ifelse(momedu=="No education", 1, 0)]
dt = dt[, edu_primary := ifelse(momedu=="Primary (1-5y)", 1, 0)]
dt = dt[, edu_secondary := ifelse(momedu=="Secondary (>5y)", 1, 0)]

row_momedu_none = mean_se_strat(
  Y_0 = dt$edu_none[dt$floor == 0],
  Y_1 = dt$edu_none[dt$floor == 1],
  Yname = "No education",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

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
# HFIAS category (binary secure/insecure)
#----------------------------------------
assert_that(names(table(is.na(dt$hfiacat)))=="FALSE") 
dt = dt[, foodsecure := ifelse(hfiacat == "Food Secure", 1, 
                               ifelse(hfiacat == "Food Insecure", 0, 999))]

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
  Y_0 = dt$elec[dt$floor == 0],
  Y_1 = dt$elec[dt$floor == 1],
  Yname = "Has electricity",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# has improved wall materials
#----------------------------------------

row_walls = mean_se_strat(
  Y_0 = dt$walls[dt$floor == 0],
  Y_1 = dt$walls[dt$floor == 1],
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
row_wardrobe = mean_se_strat(
  Y_0 = dt$asset_wardrobe[dt$floor == 0],
  Y_1 = dt$asset_wardrobe[dt$floor == 1],
  Yname = "Owns >=1 wardrobe",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_table = mean_se_strat(
  Y_0 = dt$asset_table[dt$floor == 0],
  Y_1 = dt$asset_table[dt$floor == 1],
  Yname = "Owns >=1 table",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_chair = mean_se_strat(
  Y_0 = dt$asset_chair[dt$floor == 0],
  Y_1 = dt$asset_chair[dt$floor == 1],
  Yname = "Owns >=1 chair",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_khat = mean_se_strat(
  Y_0 = dt$asset_khat[dt$floor == 0],
  Y_1 = dt$asset_khat[dt$floor == 1],
  Yname = "Owns >=1 khat",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_chouki = mean_se_strat(
  Y_0 = dt$asset_chouki[dt$floor == 0],
  Y_1 = dt$asset_chouki[dt$floor == 1],
  Yname = "Owns >=1 chouki",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_tv = mean_se_strat(
  Y_0 = dt$asset_tv[dt$floor == 0],
  Y_1 = dt$asset_tv[dt$floor == 1],
  Yname = "Owns >=1 tv",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_refrigerator = mean_se_strat(
  Y_0 = dt$asset_refrig[dt$floor == 0],
  Y_1 = dt$asset_refrig[dt$floor == 1],
  Yname = "Owns >=1 fridge",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_bike = mean_se_strat(
  Y_0 = dt$asset_bike[dt$floor == 0],
  Y_1 = dt$asset_bike[dt$floor == 1],
  Yname = "Owns >=1 bicycle",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_moto = mean_se_strat(
  Y_0 = dt$asset_moto[dt$floor == 0],
  Y_1 = dt$asset_moto[dt$floor == 1],
  Yname = "Owns >=1 motorcycle",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_sewmach = mean_se_strat(
  Y_0 = dt$asset_sewmach[dt$floor == 0],
  Y_1 = dt$asset_sewmach[dt$floor == 1],
  Yname = "Owns >=1 sewing machine",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

row_mobile = mean_se_strat(
  Y_0 = dt$asset_mobile[dt$floor == 0],
  Y_1 = dt$asset_mobile[dt$floor == 1],
  Yname = "Owns >=1 mobile phone",
  id_0 = dt$hhid[dt$floor == 0],
  id_1 = dt$hhid[dt$floor == 1]
)

#----------------------------------------
# intervention arm
#----------------------------------------
assert_that(names(table(is.na(dt$tr)))=="FALSE") 
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
# child agecat
#----------------------------------------
assert_that(names(table(is.na(dt$agecat)))=="FALSE") # Should output TRUE
dt = dt[, age_0to5 := ifelse(agecat=="0-5 years", 1, 0)]
dt = dt[, age_6to12 := ifelse(agecat =="6-12 years", 1, 0)]

row_child0to5 = mean_se_strat(
  Y_0 = dt$age_0to5[dt$floor == 0],
  Y_1 = dt$age_0to5[dt$floor == 1],
  Yname = "Child's age 0-5 years",
  id_0 = dt$personid[dt$floor == 0],
  id_1 = dt$personid[dt$floor == 1]
)

# Note KE data has 6-15 years category
row_child6to12 = mean_se_strat(
  Y_0 = dt$age_6to12[dt$floor == 0],
  Y_1 = dt$age_6to12[dt$floor == 1],
  Yname = "Child's age 6-12 years",
  id_0 = dt$personid[dt$floor == 0],
  id_1 = dt$personid[dt$floor == 1]
)

#----------------------------------------
# child mean age, years
#----------------------------------------

row_childageyr = mean_se_strat(
  Y_0 = dt$agey[dt$floor == 0],
  Y_1 = dt$agey[dt$floor == 1],
  Yname = "Child age, years",
  id_0 = dt$personid[dt$floor == 0],
  id_1 = dt$personid[dt$floor == 1]
)

#----------------------------------------
# month -- move this to table 2
#----------------------------------------

dt$svymonth <- substr(dt$svydate, 3, 5)

#----------------------------------------
# child sex
#----------------------------------------
assert_that(names(table(is.na(dt$sex)))=="FALSE") 
dt = dt[, sexmale := ifelse(sex=="male", 1, 0)]

row_childmale = mean_se_strat(
  Y_0 = dt$sexmale[dt$floor == 0],
  Y_1 = dt$sexmale[dt$floor == 1],
  Yname = "Male, %",
  id_0 = dt$personid[dt$floor == 0],
  id_1 = dt$personid[dt$floor == 1]
)


########################################################
# Make tables
########################################################

# Table of household characteristics
table = bind_rows(
  row_momage, row_momheight, row_momedu_primary, row_momedu_secondary, 
  row_Nind18, row_Nind,
  row_foodsecure, 
  row_electricity, row_walls, row_roof, row_tv, 
  row_bike, row_moto, row_mobile,
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
                  "At least some primary education", "At least some secondary education",  
  "# individuals living in compound <=18 yrs", "Total individuals living in compound",
  "Food secure", 
  "Has electricity", "Has improved wall materials", "Has improved roof material", "Owns >=1 tv", 
  "Owns >=1 bicycle", "Owns >=1 motorcycle", "Owns >=1 mobile phone",
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
compound_header = c("Compound", rep("", 4))
household_header = c("Household", rep("", 4))
# arm_header = c("Intervention assignment", rep("", 4))
child_header = c("Child", rep("", 4))

names(maternal_header) = colnames(table_wide_all_ordered)
names(compound_header) = colnames(table_wide_all_ordered)
names(household_header) = colnames(table_wide_all_ordered)
# names(arm_header) = colnames(table_wide_all_ordered)
names(child_header) = colnames(table_wide_all_ordered)

table_wide_out = bind_rows(
  maternal_header, table_wide_all_ordered[1:4,],
  compound_header, table_wide_all_ordered[5:6,],
  household_header, table_wide_all_ordered[7:14,],
  child_header, table_wide_all_ordered[15:16,]
)

table_bd_full <- table_wide_out

########################################################
# Save tables
########################################################

write.csv(table_bd_full, file=paste0(tab_path, "/table-characteristics-pos-bd.txt"), row.names=FALSE)
write.csv(table_bd_full, file=paste0(tab_path, "/table-characteristics-pos-bd.csv"), row.names=FALSE)
