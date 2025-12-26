# Clear console, environment, and garbage collect
cat("\014")       
rm(list=ls())     
gc()              

# Load required packages
library(dplyr)   
library(readr)    

# Set working directory
setwd("C:/Users/qianh/Desktop/R/cam.pla.div")
getwd()          

# Import datasets
checklist <- read_csv("./data/Raw_data/Campus_checklist_wcvp.csv")          # Main species checklist
native_directory <- read_csv("./data/Raw_data/Native_directory_wcvp.csv")   # Native species list
alien_directory <- read_csv("./data/Raw_data/Alien_directory_wcvp.csv")     # Alien species list
invasive_directory <- read_csv("./data/Raw_data/Invasive_directory_wcvp.csv") # Invasive species list

# Deduplicate checklist - keep only first record per university/species combination
checklist <- checklist %>%
  group_by(univ.links.uni.abbrev02, taxon_name) %>%  # Group by university and species
  slice(1) %>%                                       # Keep first row in each group
  ungroup()        

# Data cleaning - check taxon ranks
table(checklist$taxon_rank, useNA = "ifany")  # Show distribution of taxonomic ranks

# Filter invasive species to only include Species rank
table(invasive_directory$taxon_rank, useNA = "ifany")
invasive_directory <- invasive_directory %>% 
  filter(taxon_rank == "Species")  # Keep only species-level records

# Filter native species to only include Species rank
table(native_directory$taxon_rank, useNA = "ifany")
native_directory <- native_directory %>% 
  filter(taxon_rank == "Species")  # Keep only species-level records

# Filter alien species to only include Species rank
table(alien_directory$taxon_rank, useNA = "ifany")
alien_directory <- alien_directory %>% 
  filter(taxon_rank == "Species")  # Keep only species-level records

# Assign status based on species lists
checklist$status <- NA  # Initialize status column with NAs
checklist$status[checklist$taxon_name %in% alien_directory$taxon_name] <- "alien"
checklist$status[checklist$taxon_name %in% native_directory$taxon_name] <- "native"
checklist$status[checklist$taxon_name %in% invasive_directory$taxon_name] <- "invasive"

# Check status distribution and remove NA status records
table(checklist$status, useNA = "ifany")
checklist <- checklist %>%
  filter(!is.na(checklist$status))  # Remove rows with NA status

# Create binary flags for native and invasive status
checklist$native <- ifelse(checklist$status == "native", 1, 0)     # 1 for native, 0 otherwise
checklist$invasive <- ifelse(checklist$status == "invasive", 1, 0) # 1 for invasive, 0 otherwise
checklist$non_invasive <- ifelse(checklist$status == "alien", 1, 0) # 1 for invasive, 0 otherwise

# Remove universities 'fjsfu', 'gxsfu', 'znlykju', and 'zust' from the checklist dataset
# Fujian Normal University, 
# Guangxi Normal University, 
# Central South University of Forestry and Technology, 
# and Zhejiang University of Science and Technology, do not have invasive plants
checklist <- checklist %>%
  filter(!univ.links.uni.abbrev02 %in% c("fjsfu", "gxsfu", "znlykju", "zust"))

# Export final dataset
write.csv(checklist, file = "./data/Campus_checklist_status.csv", row.names = FALSE)