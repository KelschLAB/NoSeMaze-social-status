# --------------------------------------------------------------
# Script: combine_sequence_files_with_ELO_ratings
# Author: Jonathan Reinwald
# Date: 17.03.2025
# Description:
# This script merges dynamic ELO ratings (winner/loser ID, pre/post ELO ratings)
# with tube competition sequence data (entry/exit times, duration, etc.) 
# for each tube competition across all cohorts. # The final dataset is 
# structured for LME modeling, assessing win/loss outcomes # based on entry time 
# and ELO rating differences.
# !!! Be sure to run "myRootPath/NoSeMaze_Experiment_Social_Status/analysis/scripts/cross_cohort_summary/main_calculate_tube_metrics_R.R 
# before to calculate the ELO ratings and save the cohort_characteristic-file before !!!

# Steps:
# 1. Load required libraries and define file paths.
# 2. Extract cohort names and iterate through each cohort.
# 3. Process dynamic ELO ratings from cohort characteristics.
# 4. Process sequence files for tube competitions.
# 5. Merge datasets to create the final structured dataset.
# 6. Output summary of the merged dataset.
# --------------------------------------------------------------

# clearing
rm(list = ls())

# Load required libraries
library(dplyr)
library(readr)
library(stringr)

# Define base directory for cohort sequence files
base_dir <- "myRootPath/NoSeMaze_Experiment_Social_Status/data/processed" # --> CHANGE HERE TO YOU OWN ROOT PATH FOLDER 

# Define path to cohort characteristics file
cohort_characteristics_file <- file.path(base_dir, "cross_cohort_files", "tube", "CohortCharacteristics_TubeCompetitions_D1_End_Koptimal.Rdata")

# Check if the file exists before proceeding
if (!file.exists(cohort_characteristics_file)) {
  stop("Error: Cohort characteristics file not found!")
}

# Load cohort characteristics data
load(cohort_characteristics_file)

# Extract Koptimal/K100 from filename
cohort_type <- ifelse(grepl("Koptimal", cohort_characteristics_file), "Koptimal", "K100")

# Initialize empty lists to store logtable (ELO) and cohort data (sequence information)
logtable_list <- list()
cohort_data_list <- list()

# Extract cohort names dynamically from cohort_characteristics
# Note that cohorts 16, 20, 21 are not considered in this analysis (low number of tube competitions)
distinct_cohorts <- cohort_characteristics$cohort  # Extract unique cohort names

# Iterate through each cohort to process data
for (cohort_name in distinct_cohorts) {
  
  # Identify the cohort index for reference
  cohort_index <- which(cohort_characteristics$cohort == cohort_name)
  
  ##############################################################################
  # Processing Dynamic ELO Ratings
  
  # The logtable file is the file containing all information on the dynamic elo ratings
  # It is read from the cohort_characteristics and contains the winner ID, loser ID, 
  # pre/post ELO ratings for both of them, and the number of the event.
  
  if (!is.null(cohort_characteristics$res[[cohort_index]]$logtable)) {
    logtable <- cohort_characteristics$res[[cohort_index]]$logtable %>%
      mutate(cohort = cohort_name) %>%
      arrange(winner, loser, Date) %>%  # Ensure sorted order
      mutate(row_id = row_number())  # Assign unique row identifier
    
    logtable_list[[cohort_name]] <- logtable
  } else {
    message(paste("Skipping", cohort_name, "- No logtable found"))
  }

  ##############################################################################
  # Processing Tube Competition Sequence Files
  # The sequence files contain information on entry and exit times and durations.
  
  # Define path to the tube competitions sequence files
  cohort_path <- file.path(base_dir, cohort_name, "tube", "sequence_files", "tube_competitions")
  # Identify the relevant CSV sequence file(s)
  csv_files <- list.files(cohort_path, pattern = "^SequenceList_.*\\.csv$", full.names = TRUE)
  
  if (length(csv_files) > 0) {
    file_path <- csv_files[1]  # Select the first matching file
    cohort_data <- read_csv(file_path, show_col_types = FALSE)
    
    # Validate the presence of required columns before proceeding
    required_columns <- c("winners", "losers", "winner_entry", "loser_entry", "winner_duration", "loser_duration")
    
    if (all(required_columns %in% colnames(cohort_data))) {
      cohort_data <- cohort_data %>%
        rename(winner = winners, loser = losers) %>%  # Standardize column names
        mutate(cohort = cohort_name) %>%
        arrange(winner, loser) %>%  # Ensure sorting consistency
        mutate(row_id = row_number())  # Assign row identifier
      
      cohort_data_list[[cohort_name]] <- cohort_data
    } else {
      message(paste("Skipping", cohort_name, "- Missing required columns"))
    }
  } else {
    message(paste("Skipping", cohort_name, "- No sequence file found"))
  }
}

################################################################################
# Combine Processed Data

# Merge all logtable data
elo_data <- bind_rows(logtable_list)

# Merge all cohort sequence data
sequence_data <- bind_rows(cohort_data_list)

# Final merge: Match logtable with sequence data by cohort, winner, loser, and row_id
merged_tube_data <- elo_data %>%
  left_join(sequence_data, by = c("winner", "loser", "cohort", "row_id"))

# Output summary
summary(merged_tube_data)

################################################################################
# Saving the Processed Dataset

# Define save directory
save_dir <- file.path(base_dir, "cross_cohort_files", "tube", "merged_elo_sequence_data")

# Create directory if it doesn't exist
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# Define output filename
output_filename <- paste0("Merged_Tube_Competitions_", cohort_type, ".csv")

# Save dataset
write_csv(merged_tube_data, file.path(save_dir, output_filename))

# Print success message
message("Dataset successfully saved to: ", file.path(save_dir, output_filename))

