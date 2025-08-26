# =============================================================
# 📌 Script: calculate_adjusted_ELO_ratings.R
# Author: Jonathan Reinwald
# Date: 17.03.2025
# -------------------------------------------------------------
# This script calculates dynamic ELO ratings for tube competitions
# across multiple cohorts. It includes both:
# - Standard Elo Rating (uncorrected)
# - Corrected Elo Rating (adjusted for entry time differences)
#
# Uses a flexible K-factor:
# - "Koptimal": Uses cohort-specific optimal K values.
# - "K100": Uses a fixed K-value of 100 for all cohorts.
# =============================================================

# 📌 Clear environment before execution
rm(list = ls())

# =============================================================
# 🚀 1. Define Parameters & Load Libraries
# =============================================================

# ✅ Define the ELO rating method (choose "Koptimal" or "K100")
cohort_type <- "Koptimal" #optimal"  # Change to "K100" if using fixed K=100

# ✅ Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(aniDom)
library(EloRating)

# ✅ Define the base directory for data
main_dir <- "myRootPath/NoSeMaze_Experiment/" # --> CHANGE HERE TO YOU OWN ROOT PATH FOLDER 
base_dir <- file.path(main_dir,"data", "processed")

# ✅ Define path to cohort characteristics file (contains optimal K values)
cohort_characteristics_filename <- paste0("CohortCharacteristics_TubeCompetitions_D1_End_", cohort_type, ".Rdata")
cohort_characteristics_file <- file.path(base_dir, "cross_cohort_files", "tube", cohort_characteristics_filename)

# ✅ Ensure cohort characteristics file exists
if (!file.exists(cohort_characteristics_file)) {
  stop("Error: Cohort characteristics file not found!")
}

# ✅ Load cohort characteristics
load(cohort_characteristics_file)

# ✅ Define cohorts (excluding cohort 16, 20, 21 due to low competition numbers)
cohort_names <- c("cohort01", "cohort02", "cohort03", "cohort04", "cohort05",
                  "cohort06", "cohort07", "cohort08", "cohort09", "cohort10",
                  "cohort11", "cohort12", "cohort13", "cohort14", "cohort15",
                  "cohort17", "cohort18", "cohort19") 

# ✅ Define initial ELO rating
initial_ELO <- 0

# =============================================================
# 📂 2. Load Scaling Parameters
# =============================================================

# ✅ Define path to scaling parameters file
data_dir <- file.path(base_dir, "cross_cohort_files", "tube", "merged_elo_sequence_data")
input_filename <- paste0("scaling_parameters_", cohort_type, ".csv")
file_path <- file.path(data_dir, input_filename)

# ✅ Load scaling parameters
if (file.exists(file_path)) {
  scaling_results <- read_csv(file_path, show_col_types = FALSE)
} else {
  stop("Error: Scaling parameters file not found!")
}


# =============================================================
# 🔄 3. Loop Through Each Cohort
# =============================================================

# ✅ Initialize list to store ELO results
all_ELO_results <- list()

# Loop through each cohort
for (cohort_name in cohort_names) {

    # ✅ Define path to tube competition sequence files
  cohort_path <- file.path(base_dir, cohort_name, "tube", "sequence_files", "tube_competitions")
  
  # ✅ Find and load CSV files for each cohort
  csv_files <- list.files(cohort_path, pattern = "^SequenceList_.*\\.csv$", full.names = TRUE)
  
  # Load all CSV files found
  cohort_data_list <- list()
  for (file_path in csv_files) {
    cohort_data <- read_csv(file_path, show_col_types = FALSE) %>%
      mutate(cohort = cohort_name)  # Add cohort identifier
    cohort_data_list[[length(cohort_data_list) + 1]] <- cohort_data
  }
  
  # ✅ Assign K-factor (cohort-specific or fixed)
  if (cohort_type == "Koptimal") {
    K_factor <- cohort_characteristics$optimal_K[cohort_characteristics$cohort == cohort_name]
  } else if (cohort_type == "K100") {
    K_factor <- 100
  }
  
  # ✅ Skip cohort if no data was found
  if (length(cohort_data_list) == 0) next  
  
  # ✅ Combine all CSV data for the cohort
  sequence_data <- bind_rows(cohort_data_list)
  
  # ✅ Calculate entry time difference
  sequence_data <- sequence_data %>%
    mutate(day = as.Date(day, format="%d.%m.%Y"),  
           entry_diff = winner_entry - loser_entry)  # Calculate entry time difference
  
  # ✅ Get unique players in the cohort
  players <- unique(c(sequence_data$winners, sequence_data$losers))
  
  # ✅ Initialize ELO ratings
  ELO_ratings <- data.frame(Player_ID = players, 
                            cohort = cohort_name,
                            ELO_uncorrected = rep(initial_ELO, length(players)),
                            ELO_corrected = rep(initial_ELO, length(players)))
  
  # ✅ Extract winners, losers, and dates for ELO calculation
  winners <- sequence_data$winners
  losers <- sequence_data$losers
  all_dates <- as.Date(sequence_data$day)
  
  # ✅ Compute ELO ratings using the EloRating package
  res <- elo.seq(winners, losers, all_dates, draw = NULL, presence = NULL, startvalue = 0,
                 k = K_factor, normprob = TRUE, init = "average", intensity = NULL,
                 iterate = 0, runcheck = TRUE, progressbar = FALSE)

  # ✅ Extract the final ELO scores at the last recorded date
  # 1. EloRating Package
  traditional_elo_1 <- extract_elo(res, extractdate = tail(all_dates, 1))
  # 2. AniDom Package
  non_randomized_ELOscore <- elo_scores(winners, losers, identities = NULL, sigmoid.param = 1/100,
                                            K = K_factor, init.score = 0, randomise = FALSE, n.rands = 1000,
                                            return.as.ranks = FALSE, return.trajectories = FALSE, dates = NULL)
  # Store row names first
  non_randomized_names <- rownames(non_randomized_ELOscore)
    # Convert the matrix column to a numeric vector
  traditional_elo_2 <- as.numeric(non_randomized_ELOscore[,1])  
  # Reassign the stored names
  names(traditional_elo_2) <- non_randomized_names  
  # Ensure both have the same names and order them correctly
  common_ids <- intersect(names(traditional_elo_1), names(traditional_elo_2))
    # Subset and reorder both vectors to match the same order
  elo1 <- traditional_elo_1[common_ids]
  elo2 <- traditional_elo_2[common_ids]
  
  # ✅ Convert traditional Elo results to a DataFrame and filter only relevant players
  traditional_elo_1_df <- data.frame(
    Player_ID = names(traditional_elo_1),
    ELO_traditional = as.numeric(traditional_elo_1)
  ) %>%
    filter(Player_ID %in% players)  # **Ensure only current cohort players are included**
  traditional_elo_2_df <- data.frame(
    Player_ID = names(traditional_elo_2),
    ELO_traditional = as.numeric(traditional_elo_2)
  ) %>%
    filter(Player_ID %in% players)  # **Ensure only current cohort players are included**
  
  # =============================================================
  # 📌 Extract Scaling Parameters for Corrections
  # =============================================================
  
  # ✅ Extract relevant scaling parameters from the dataset
  max_entry_diff <- scaling_results %>% filter(Parameter == "max_entry_diff_threshold") %>% pull(Value)
  max_duration_threshold <- scaling_results %>% filter(Parameter == "max_duration_threshold") %>% pull(Value)
  mean_ELO_diff <- scaling_results %>% filter(Parameter == "mean_ELO_diff") %>% pull(Value)
  sd_ELO_diff <- scaling_results %>% filter(Parameter == "sd_ELO_diff") %>% pull(Value)
  mean_entry_diff <- scaling_results %>% filter(Parameter == "mean_entry_diff") %>% pull(Value)
  sd_entry_diff <- scaling_results %>% filter(Parameter == "sd_entry_diff") %>% pull(Value)
  entry_diff_scaled_estimate <- scaling_results %>% filter(Parameter == "entry_diff_scaled_estimate") %>% pull(Value)
  
  # ✅ Compute scaled K-factor based on ELO standard deviation
  K_factor_scaled <- K_factor / sd_ELO_diff  
  
  # =============================================================
  # 📌 Define Function for Corrected Win Probability
  # =============================================================
  
  calculate_win_probability_corrected <- function(ELO_winner, ELO_loser, entry_diff) {
    # Standardize ELO difference
    ELO_diff_scaled <- (ELO_winner - ELO_loser - mean_ELO_diff) / sd_ELO_diff
    expected_win <- 1 / (1 + 10^(-ELO_diff_scaled / 400))
    
    # Compute correction based on entry time difference
    beta_entry <- entry_diff_scaled_estimate
    entry_diff_scaled <- (entry_diff - mean_entry_diff) / sd_entry_diff  
    adjusted_expected_win <- expected_win + (beta_entry * entry_diff_scaled)
    
    # Ensure probability remains between 0.01 and 0.99
    return(max(min(adjusted_expected_win, 0.99), 0.01))
  }
  
  # =============================================================
  # 🔄 Process Competitions for This Cohort
  # =============================================================
  
  # Process Competitions for this cohort
  for (match in 1:nrow(sequence_data)) {
    
    winner <- sequence_data$winners[match]
    loser <- sequence_data$losers[match]
    entry_diff <- sequence_data$entry_diff[match]
    
    # ✅ Cap entry_diff within allowed thresholds
    if (entry_diff > max_entry_diff) entry_diff <- max_entry_diff
    if (entry_diff < -max_entry_diff) entry_diff <- -max_entry_diff
    
    # ✅ Retrieve ELO ratings before match
    ELO_winner_uncorrected <- ELO_ratings$ELO_uncorrected[ELO_ratings$Player_ID == winner]
    ELO_loser_uncorrected <- ELO_ratings$ELO_uncorrected[ELO_ratings$Player_ID == loser]
    ELO_winner_corrected <- ELO_ratings$ELO_corrected[ELO_ratings$Player_ID == winner]
    ELO_loser_corrected <- ELO_ratings$ELO_corrected[ELO_ratings$Player_ID == loser]
    
    # ✅ Compute expected win probabilities (uncorrected and corrected)
    expected_win_uncorrected <- 1 / (1 + 10^(-(ELO_winner_uncorrected - ELO_loser_uncorrected) / 400))
    expected_win_corrected <- calculate_win_probability_corrected(ELO_winner_corrected, ELO_loser_corrected, entry_diff)
    
    # ✅ Update Uncorrected ELO Ratings
    ELO_ratings$ELO_uncorrected[ELO_ratings$Player_ID == winner] <- 
      ELO_winner_uncorrected + K_factor * (1 - expected_win_uncorrected)
    
    ELO_ratings$ELO_uncorrected[ELO_ratings$Player_ID == loser] <- 
      ELO_loser_uncorrected + K_factor * (0 - (1 - expected_win_uncorrected))
    
    # ✅ Update Corrected ELO Ratings
    ELO_ratings$ELO_corrected[ELO_ratings$Player_ID == winner] <-
      ELO_winner_corrected + K_factor_scaled * (1 - expected_win_corrected)

    ELO_ratings$ELO_corrected[ELO_ratings$Player_ID == loser] <-
      ELO_loser_corrected + K_factor_scaled * (0 - (1 - expected_win_corrected))
    
    # scaling_factor <- sd_ELO_diff  # Undo the original downscaling
    # 
    # # ✅ Update Corrected ELO Ratings (Properly Scaled)
    # ELO_ratings$ELO_corrected[ELO_ratings$Player_ID == winner] <- 
    #   ELO_winner_corrected + (K_factor_scaled  * scaling_factor) * (1 - expected_win_corrected)
    # 
    # ELO_ratings$ELO_corrected[ELO_ratings$Player_ID == loser] <- 
    #   ELO_loser_corrected + (K_factor_scaled * scaling_factor) * (0 - (1 - expected_win_corrected))
    
  }
  
  # =============================================================
  # 📌 Store Results for This Cohort
  # =============================================================
  
  # Merge the two traditional ELO score dataframes
  traditional_elo_combined <- full_join(traditional_elo_1_df, traditional_elo_2_df, by = "Player_ID")
  
  # Rename columns explicitly to avoid confusion
  colnames(traditional_elo_combined) <- c("Player_ID", "Traditional_ELO_1", "Traditional_ELO_2")
  
  # Rescale corrected ELOs to have the same mean and SD as uncorrected ELOs
  unc_mean <- mean(traditional_elo_combined$Traditional_ELO_2, na.rm = TRUE)
  unc_sd <- sd(traditional_elo_combined$Traditional_ELO_2, na.rm = TRUE)
  cor_mean <- mean(ELO_ratings$ELO_corrected, na.rm = TRUE)
  cor_sd <- sd(ELO_ratings$ELO_corrected, na.rm = TRUE)
  
  # Apply linear rescaling
  ELO_ratings$ELO_corrected_rescaled <- ((ELO_ratings$ELO_corrected - cor_mean) / cor_sd) * unc_sd + unc_mean

  # ✅ Perge the combined ELO scores with ELO_ratings
  all_ELO_results[[cohort_name]] <- ELO_ratings %>%
    left_join(traditional_elo_combined, by = "Player_ID")
  
  # ✅ Print results for debugging (optional)
  print(all_ELO_results[[cohort_name]])
}

# ✅ Define results directory
save_dir <- file.path(base_dir, "cross_cohort_files", "tube", "merged_elo_sequence_data")


# ✅ Create directory if it doesn’t exist
if (!dir.exists(save_dir)) {
  dir.create(save_dir, recursive = TRUE)
}

# ✅ Define Save Path with cohort_type in filename
save_filename <- paste0("Final_ELO_ratings_by_cohort_", cohort_type, ".csv")
save_path <- file.path(save_dir, save_filename)

# ✅ Save as CSV
final_ELO_ratings <- bind_rows(all_ELO_results)
write.csv(final_ELO_ratings, save_path, row.names = FALSE)

# 📊 Check correlation
elo_correlation <- cor(final_ELO_ratings$ELO_uncorrected, final_ELO_ratings$ELO_corrected, use = "complete.obs")
print(paste("Correlation between corrected and uncorrected ELO:", elo_correlation))

# ✅ Define results directory
results_dir <- file.path(main_dir, "results", "figures", "cross_cohort", "adjusted_for_entry_difference")


# ✅ Create directory if it doesn’t exist
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Update correlation
elo_correlation <- cor(final_ELO_ratings$ELO_uncorrected, final_ELO_ratings$ELO_corrected_rescaled, use = "complete.obs")

# Compute symmetric max absolute value from both ELO types
elo_range <- range(
  final_ELO_ratings$ELO_uncorrected,
  final_ELO_ratings$ELO_corrected_rescaled,
  na.rm = TRUE
)
elo_abs_max <- max(abs(elo_range))
elo_lim <- ceiling(elo_abs_max / 50) * 50  # Rounded axis limit

# Create plot
elo_plot <- ggplot(final_ELO_ratings, aes(x = Traditional_ELO_2, y = ELO_corrected_rescaled)) +
  geom_point(alpha = 0.6, color = "blue") +
  geom_abline(slope = 1, intercept = 0, color = "gray40", linetype = "dotted") +
  geom_smooth(method = "lm", color = "red", linetype = "dashed") +
  coord_fixed(ratio = 1, xlim = c(-elo_lim, elo_lim), ylim = c(-elo_lim, elo_lim)) +
  labs(
    title = "Corrected vs. Uncorrected ELO Ratings (Rescaled)",
    subtitle = paste("Pearson r =", round(elo_correlation, 3)),
    x = "Uncorrected ELO",
    y = "Corrected ELO (Rescaled)"
  ) +
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14)
  )



# ✅ Define file path for saving the figure
plot_name <- paste0("ELO_Correlation_Corrected_vs_Uncorrected_", cohort_type, ".pdf")
plot_filename <- file.path(results_dir, plot_name)

# ✅ Save plot as a high-resolution PNG
ggsave(plot_filename, elo_plot, width = 8, height = 6, dpi = 300)

# ✅ Print confirmation message
message("ELO correlation plot saved to: ", plot_filename)

# 📁 Save source data used for ELO correlation analysis

# Create a clean data frame of the two columns used in the correlation
elo_correlation_data <- final_ELO_ratings %>%
  select(Player_ID, cohort, ELO_uncorrected, ELO_corrected_rescaled) %>%
  filter(!is.na(ELO_uncorrected) & !is.na(ELO_corrected_rescaled))  # Ensure complete cases only

# Define save path
source_data_path <- file.path(results_dir, paste0("SourceData_ELOCorrelation_", cohort_type, ".csv"))

# Write CSV
write.csv(elo_correlation_data, source_data_path, row.names = FALSE)

# Confirm
message("✅ Source data for ELO correlation saved to: ", source_data_path)
