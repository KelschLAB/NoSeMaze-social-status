# Main script to analyze tube data from NoSeMaze cohorts
# Author: Jonathan Reinwald, Revised: 17.06.2024
# Be sure to use a relatively recent R-version, e.g., /zi/apps/bin/R_4.2.3

# Pre-Clearing
rm(list = ls())

# Libraries
library(dplyr)
library(readr)
# in the subscripts: library(aniDom)
# in the subscripts: library(EloRating)

# Predefine K-value and day ranges
day_ranges <- c("D1_End", "D1_14", "D1_21", "Last14", "D1_7", "D8_14", "D15_21")
# day_ranges <- c("D1_7", "D8_14", "D15_21", "Last14")
# K <- 100
K <- 'optimal'

# Specify directories
base_dir <- "myRootPath/NoSeMaze_Experiment" # --> REPLACE WITH YOUR OWN ROOT DIRECTORY
data_dir <- file.path(base_dir, "data/processed")
config_file <- file.path(base_dir, "config/cohorts_info.csv")

# Load cohort configuration
cohorts_tbl <- read_csv(config_file)
selected_cohorts <- cohorts_tbl %>%
  filter(use_tube == 1) %>%
  pull(cohort)

# Load source script
source(file.path(base_dir, "/src/tube/cross_cohort_summary/calculate_tube_metrics_R.R"))

# Define types of competitions
data_types <- c("TubeCompetitions", "TubeChasings") #

# Ensure directories exist or create them
output_dirs <- list(
  cohort_files_dir = file.path(data_dir, "cross_cohort_files/tube"),
  correlation_results_dir = file.path(base_dir, "results/figures/cross_cohort/stability_across_metrics/tube")
)


# Ensure directories exist or create them
lapply(output_dirs, function(dir) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
})

# Loop over competition types and day ranges
for (data_type in data_types) {
  for (day_range_name in day_ranges) {

    # Ensure correct path separators
    data_dir <- normalizePath(data_dir, winslash = "/", mustWork = FALSE)
    
    # Get list of files for the current competition type
    filelist <- list.files(
      path = data_dir,
      pattern = paste0("^SequenceList_.*_", data_type, "\\.csv$"),
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE  # Make it case-insensitive for Linux compatibility
    )
    
    # Exclude cohorts by date range
    filtered_by_date <- filter_files_by_date_range(filelist, day_range_name, cohorts_tbl)

    # Exclude cohorts not in the selected cohorts list
    filtered_filelist <- Filter(function(file) {
      cohort_name <- sub('.*SequenceList_(.*?)_days.*', '\\1', basename(file))
      cohort_name %in% selected_cohorts
    }, filtered_by_date)

    # Run calculations
    results <- calculate_cohort_characteristics(filtered_filelist, day_range_name, K, cohorts_tbl)
    
    browser()
    # Direct assignment (necessary for saving)
    cohort_characteristics <- results$cohort_characteristics
    correlation_results <- results$correlation_results
    
    # Save cohort_characteristics as RData and CSV
    save(cohort_characteristics, file = file.path(output_dirs$cohort_files_dir, 
                                                  paste0("CohortCharacteristics_", data_type, "_", day_range_name, "_K", K, ".RData")))    
    cohort_characteristics <- cohort_characteristics %>%
      mutate(across(where(is.list), ~sapply(., function(x) if (length(x) == 1) x else paste(x, collapse = ";"))))
    write.csv(cohort_characteristics, file = file.path(output_dirs$cohort_files_dir, 
                                                       paste0("CohortCharacteristics_", data_type, "_", day_range_name, "_K", K, ".csv")), 
              row.names = FALSE)    
    
    # Save correlation_results as RData and CSV
    save(correlation_results, file = file.path(output_dirs$correlation_results_dir, 
                                               paste0("CorrelationResults_", data_type, "_", day_range_name, "_K", K, ".RData")))
    correlation_results <- correlation_results %>%
      mutate(across(where(is.list), ~sapply(., function(x) if (length(x) == 1) x else paste(x, collapse = ";"))))
    write.csv(correlation_results, file = file.path(output_dirs$correlation_results_dir, 
                                                    paste0("CorrelationResults_", data_type, "_", day_range_name, "_K", K, ".csv")), 
              row.names = FALSE)
    
  }
}