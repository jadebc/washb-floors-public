#######################################
# WASH Benefits STH finished floor analysis

# Objective: this code generates a full table 1 - joined BD and KE
#######################################


# Configure directories, load libraries and base functions
source(paste0(here::here(), "/0-config.R"))
library(dplyr)
library(data.table)

# Load in Bangladesh table 
bdtable <- read.table(file=paste0(tab_path, "/table-characteristics-bd.txt"), sep = ",", header = TRUE)

# Load in Kenya table
ketable <- read.table(file=paste0(tab_path, "/table-characteristics-ke.txt"), sep=",", header = TRUE)

# Merge tables together on variable
fulltable1 <- merge(bdtable, ketable, 
                    by.x="Variable", by.y="Variable", 
                    all.y=TRUE, all.x=TRUE)

# rename columns
colnames(fulltable1) = c(
  "Variable",
  "(BD) N, finished", "(BD) Results, finished",
  "(BD) N, unfinished", "(BD) Results, unfinished",
  "(KE) N, finished", "(KE) Results, finished",
  "(KE) N, unfinished", "(KE)Results, unfinished"
)

# Resort by variable name
full_outcome_list <- c("Maternal",
                       "Mother's age, years", "Mother's height, cm", 
                       "At least some primary education", "At least some secondary education",
                       "Compound",
                       "# individuals living in compound <=18 yrs", "Total individuals living in compound",
                       "Household", "Food secure", 
                       "Has electricity", "Has improved wall materials", "Has improved roof material", 
                       "Owns >=1 tv", "Owns >=1 bicycle", "Owns >=1 motorcycle", "Owns >=1 mobile phone", 
                       "Child", "Child age, years", "Male, %"
)

fulltable1_ordered <- fulltable1 %>%
  arrange(match(Variable, full_outcome_list))

# Clear cells with NA, change to blanks
fulltable1_ordered <- sapply(fulltable1_ordered, as.character)
fulltable1_ordered[is.na(fulltable1_ordered)] <- ""

# Eliminate duplicates on merge
fulltable1_ordered <- fulltable1_ordered[!duplicated(fulltable1_ordered),]

########################################################
# Save tables
########################################################

write.csv(fulltable1_ordered, file=paste0(tab_path, "/table-characteristics-full.csv"), row.names=FALSE)
