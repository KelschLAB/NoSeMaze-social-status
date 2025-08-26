# Script to boxplot stability measures and correlation results for tube data
# Author: Jonathan Reinwald, Revised: 16.01.2025

# Pre-clearing
rm(list = ls())

# Specify directories
base_dir <- "myRootPath/NoSeMaze_Experiment" # --> replace with your own root directory
config_file <- file.path(base_dir, "/config/cohorts_info.csv")

# Source script
source(file.path(base_dir, "/src/analysis/stability_measures_tube_metrics/plot_stability_measures_tube_R.R"))

# Define configurations
datasets <- list(
  cohort_characteristics = list(
    variables = c("transitivity_pt", "steepness", "linearity_h2", "stabilityIndex", "numberofevents", "uncertainty_rep", "optimal_K"),
    input_data_dir = file.path(base_dir, "/data/processed/cross_cohort_files/tube"),
    file_prefix = "CohortCharacteristics",
    plot_dir = file.path(base_dir, "results/figures/cross_cohort/boxplots_cohort_characteristics/tube"),
    data_variable = "cohort_characteristics"
  ),
  correlation_results = list(
    variables = c("randELOtoELO_estimate_p", "randELOtoDS_estimate_p", "ELOtoDS_estimate_p", "ELOaniDomtoELOEloRating_estimate_p"),
    input_data_dir = file.path(base_dir, "/results/figures/cross_cohort/stability_across_metrics/tube"),
    file_prefix = "CorrelationResults",
    plot_dir = file.path(base_dir, "results/figures/cross_cohort/stability_across_metrics/tube"),
    data_variable = "correlation_results"
  )
)
day_ranges <- c("D1_End", "D1_14", "D1_21")
data_types <- c("TubeCompetitions", "TubeChasings")
script_name <- "Script: C: ~/NoSeMaze_Experiment_Social_Status/analysis/scripts/stability_measures/tube/boxplot_stability_measures_tube_metrics_R.R"
K <- 100 #'optimal'

# Load cohort configuration
cohorts_tbl <- read_csv(config_file)
selected_cohorts <- cohorts_tbl %>% filter(use_tube == 1) %>% pull(cohort)

# Custom color palette
cohort_colors <- colorRampPalette(brewer.pal(12, "Set3"))(18)

# Loop over datasets, data types, and day ranges
for (dataset_name in names(datasets)) {
  dataset_config <- datasets[[dataset_name]]
  input_data_dir <- dataset_config$input_data_dir
  variables <- dataset_config$variables
  plot_subdir <- dataset_config$plot_subdir
  file_prefix <- dataset_config$file_prefix
  data_variable <- dataset_config$data_variable
  
  # Ensure output directory exists for the dataset
  dataset_plot_dir <- dataset_config$plot_dir
  ensure_directory(dataset_plot_dir)
  
  for (data_type in data_types) {
    for (day_range_name in day_ranges) {
      
      file_path <- file.path(input_data_dir, paste0(file_prefix, "_", data_type, "_", day_range_name, ".RData"))
      
      if (file.exists(file_path)) {
        load(file_path)  # Load the data (assigns either `cohort_characteristics` or `correlation_results`)
        
        # Dynamically use the correct data variable
        data_to_plot <- get(data_variable)  # Extract the appropriate data frame (e.g., `cohort_characteristics` or `correlation_results`)
        
        # Generate and save plots
        plots <- generate_variable_boxplots(data_to_plot, variables, cohort_colors)
        save_combined_plots(
          plots,data_to_plot,
          file.path(dataset_plot_dir, paste0("Plots_", dataset_name, "_", data_type, "_", day_range_name, "_K", K)),
          script_name
        )
      } else {
        message("File not found: ", file_path)
      }
    }
  }
}
