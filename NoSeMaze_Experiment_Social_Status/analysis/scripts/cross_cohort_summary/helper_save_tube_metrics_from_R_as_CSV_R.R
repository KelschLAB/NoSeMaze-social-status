# Helper function: 

# Pre-Clearing
rm(list = ls())

# Libraries
library(dplyr)

# Directory containing the RData files
base_dir <- "myRootPath/NoSeMaze_Experiment" # --> REPLACE WITH YOUR OWN ROOT DIRECTORY
data_dir <- file.path(base_dir, "data/processed/cross_cohort_files/tube")
output_dir <- file.path(base_dir, "data/interim/cross_cohort_files/tube")
data_types <- c("TubeCompetitions", "TubeChasings")
day_ranges <- c("D1_7", "D8_14", "D15_21", "D1_14", "D1_21", "D1_End", "Last14")
K = "optimal" # or 100

# Create the data directory if it does not exist
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# Initialize an empty list to store results
all_data <- list()

# Loop over competition types and day ranges
for (data_type in data_types) {
  for (day_range_name in day_ranges) {
    
    # List all RData files in the directory for the current type and day range
    rdata_files <- list.files(
      data_dir, 
      pattern = paste0("CohortCharacteristics_", data_type, "_", day_range_name, "_K", K, ".*\\.RData"), 
      full.names = TRUE
    )
    # Check if no files are found and stop execution
    if (length(rdata_files) == 0) {
      stop(sprintf("Error: No RData files found in %s matching pattern: CohortCharacteristics_%s_%s_K%s.*.RData", 
                   data_dir, data_type, day_range_name, K))
    }
    
    # Output
    output_file <- file.path(output_dir, paste0("CohortCharacteristics_", data_type, "_", day_range_name, "_K", K, "_combined.csv"))
    
    # Iterate over all RData files
    for (file in rdata_files) {
      # Initialize an empty list to store results for the current file
      all_data <- list()
      
      # Load the RData file
      load(file)
      
      # Extract the required data
      cohort_names <- cohort_characteristics$cohort
      sorted_randELOscore <- cohort_characteristics$sorted_randELOscore
      sorted_non_randELOscore <- cohort_characteristics$sorted_non_randELOscore
      sorted_DS <- cohort_characteristics$sorted_DS
      ids_list <- lapply(sorted_DS, rownames) # Extract row names for each sorted_DS
      
      # Check if the data was extracted correctly
      if (is.null(cohort_names) || is.null(sorted_randELOscore) || is.null(sorted_non_randELOscore) || is.null(sorted_DS)) {
        warning(sprintf("File %s could not be processed: Missing data", file))
        next
      }

      # Combine ids_list, sorted_DS, sorted_randELO, and sorted_non_randELO into separate columns
      for (i in seq_along(sorted_DS)) {
        ids <- ids_list[[i]]
        ds <- as.vector(sorted_DS[[i]])
        randELO <- as.vector(sorted_randELOscore[[i]])
        non_randELO <- as.vector(sorted_non_randELOscore[[i]])
        
        n_rows <- min(length(ids), length(ds), length(randELO), length(non_randELO))
        cohort <- rep(cohort_names[[i]], n_rows)
        
          # Check if all vectors have the same length
        if (n_rows == 0) {
          warning(sprintf("Empty data in file %s, cohort %s, iteration %d. Skipping.", file, cohort, i))
          next
        }

        # Create the data frame
        cohort_data <- data.frame(
          ID = ids,
          cohort = cohort,
          sorted_DS = ds,
          sorted_randELOscore = randELO,
          sorted_non_randELOscore = non_randELO
        )
        # Add the table to the overall list
        all_data[[length(all_data) + 1]] <- cohort_data
        
      }
    }
    # Combine all tables into a single large table
    final_table <- bind_rows(all_data)
    
    # Save the table as a CSV
    write.csv(final_table, output_file, row.names = FALSE)
    
    cat("The combined table has been successfully saved to", output_file, "\n")
  }
}