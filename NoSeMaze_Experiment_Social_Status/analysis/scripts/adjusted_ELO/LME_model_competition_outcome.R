# =============================================================
# Script: LME_model_competition_outcome.R
# Author: Jonathan Reinwald
# Date: 17.03.2025
# -------------------------------------------------------------
# This script loads merged tube competition data, processes it,
# and applies a Generalized Linear Mixed Model (GLMM) to analyze
# how ELO difference and entry time influence winning probability.
# =============================================================

# 📌 Clear environment before execution
rm(list = ls())

# 📦 Load Required Libraries
library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(lme4)
library(ggplot2)
library(performance)  # For R² and model comparison

# -------------------------------------------------------------
# 🗂 1. Predefinitions 
# -------------------------------------------------------------

# Definition of K for the ELO rating
cohort_type <- "Koptimal" #"K100"

# Define flexible threshold values for outlier handling
max_entry_diff <- 10 # Threshold for maximal entry_diff
max_duration <- 15 # Threshold for durations

# -------------------------------------------------------------
# 🗂 2. Load Preprocessed Data (Merged Tube Competition Dataset)
# -------------------------------------------------------------
# Define base directory for cohort sequence files
base_dir <- "myRootPath/NoSeMaze_Experiment/data/processed"  # --> CHANGE HERE TO YOU OWN ROOT PATH FOLDER 
# Define data directory
data_dir <- file.path(base_dir, "cross_cohort_files", "tube", "merged_elo_sequence_data")
# Identify the correct file (Koptimal/K100)
input_filename <- paste0("Merged_Tube_Competitions_", cohort_type, ".csv")
file_path <- file.path(data_dir, input_filename)
# Load dataset
if (file.exists(file_path)) {
  df <- read_csv(file_path, show_col_types = FALSE)
} else {
  stop("Error: Merged dataset not found. Ensure it was generated first.")
}

# -------------------------------------------------------------
# 🔍 3. Function to Replace Outliers with NA or max_entry_diff/max_duration values
# -------------------------------------------------------------
# Function to replace values exceeding thresholds with NA and count separately
replace_outlier_values <- function(data) {
  outlier_counts <- list()
  
  # Count positive and negative outliers separately for entry_diff
  outlier_counts$entry_diff_positive <- sum(data$entry_diff > max_entry_diff, na.rm = TRUE)
  outlier_counts$entry_diff_negative <- sum(data$entry_diff < -max_entry_diff, na.rm = TRUE)
  outlier_counts$duration_exceeding <- sum(data$duration > max_duration, na.rm = TRUE)
  outlier_counts$opponent_duration_exceeding <- sum(data$opponent_duration > max_duration, na.rm = TRUE)
  
  # Replace extreme entry_diff values with NA
  # data$entry_diff[abs(data$entry_diff) > max_entry_diff] <- NA
  data$entry_diff[data$entry_diff > max_entry_diff] <- max_entry_diff
  data$entry_diff[data$entry_diff < -1*max_entry_diff] <- -1 * max_entry_diff
  
  # Replace extreme duration values with NA
  # data$duration[data$duration > max_duration] <- NA
  # data$opponent_duration[data$opponent_duration > max_duration] <- NA
  data$duration[data$duration > max_duration] <- max_duration
  data$opponent_duration[data$opponent_duration > max_duration] <- max_duration
  
  
  # Print summary of replacements
  print("Number of values replaced with NA:")
  print(outlier_counts)
  
  return(data)
}

# -------------------------------------------------------------
# 🔄 4. Reshape Data into Long Format & Compute Key Variables
# -------------------------------------------------------------
# Compute ELO_diff before converting to long format
df <- df %>%
  mutate(ELO_diff = Apre - Bpre)  # Compute ELO difference

df <- df %>%
  mutate(
    Opponent_ID_winner = loser,  # Opponent for the winner is the loser
    Opponent_ID_loser = winner   # Opponent for the loser is the winner
  )

# Convert to long format & fix variables
df_long <- df %>%
  pivot_longer(cols = c(winner, loser), names_to = "role", values_to = "Player_ID") %>%
  mutate(
    Win = ifelse(role == "winner", 1, 0),  # 1 for winner, 0 for loser
    Opponent_ID = ifelse(role == "winner", Opponent_ID_winner, Opponent_ID_loser),  # Correct opponent ID assignment
    ELO_pre = ifelse(role == "winner", Apre, Bpre),  # Individual ELO before match
    entry_time = ifelse(role == "winner", winner_entry, loser_entry),  # Individual entry time
    duration = ifelse(role == "winner", winner_duration, loser_duration),  # Individual duration
    opponent_entry_time = ifelse(role == "winner", loser_entry, winner_entry),  # Opponent's entry time
    opponent_duration = ifelse(role == "winner", loser_duration, winner_duration)  # Opponent's duration
  ) %>%
  mutate(
    # ✅ Flip ELO_diff for losers
    ELO_diff = ifelse(role == "winner", Apre - Bpre, Bpre - Apre),
    
    # ✅ Flip entry_diff so that positive means the individual entered first
    entry_diff = ifelse(role == "winner", winner_entry - loser_entry, loser_entry - winner_entry)
  ) %>%
  select(Date, Player_ID, Opponent_ID, cohort, Win, ELO_pre, ELO_diff, entry_time, entry_diff, duration, opponent_duration)

# ✅ Replace values exceeding thresholds with NA  or max values 
df_long <- replace_outlier_values(df_long)

# -------------------------------------------------------------
# 📊 5. Compute Last ELO Score for Each Player
# -------------------------------------------------------------

# Step 1: Reshape df to long format to track individual player ELOs
df_long_ELO <- df %>%
  pivot_longer(cols = c(winner, loser), names_to = "role", values_to = "Player_ID") %>%
  mutate(ELO_post = ifelse(role == "winner", Apost, Bpost))  # Get the correct post-match ELO

# Step 2: Get the last recorded ELO for each player (based on last `day`)
last_ELO_lookup <- df_long_ELO %>%
  group_by(Player_ID) %>%
  filter(day == max(day)) %>%  # Keep only the last competition day for each player
  summarise(last_ELO = ELO_post[which.max(day)], .groups = "drop")  # Take ELO from that match

# Step 3: Merge last_ELO with df_long
df_long <- df_long %>%
  left_join(last_ELO_lookup, by = "Player_ID") %>%
  mutate(last_ELO_scaled = scale(last_ELO))  # Standardize for modeling

# -------------------------------------------------------------
# 🔬 6. Standardize Predictors for Model Stability
# -------------------------------------------------------------
# ✅ Standardize predictors for stability (after replacing outliers )
df_long <- df_long %>%
  mutate(
    ELO_diff_scaled = scale(ELO_diff),
    entry_diff_scaled = scale(entry_diff),
    duration_scaled = scale(duration),
    opponent_duration_scaled = scale(opponent_duration),
    ELO_pre_scaled = scale(ELO_pre),
    entry_time_scaled = scale(entry_time)
  )

# -------------------------------------------------------------
# 📉 7. Fit Generalized Linear Mixed Model (GLMM)
# -------------------------------------------------------------
model_long <- glmer(Win ~ ELO_diff_scaled + entry_diff_scaled + (1 | Player_ID),
                    data = df_long, family = binomial,
                    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)))

# Print model summary
summary(model_long)

# -------------------------------------------------------------
# 🏆 8. Evaluate Model Performance and Save Important Results
# -------------------------------------------------------------

# Compute R² and AIC for model comparison
r2(model_long)
AIC(model_long)

# 📌 Extract Fixed Effect Estimate for entry_diff_scaled
entry_diff_scaled_estimate <- fixef(model_long)["entry_diff_scaled"]

# 📌 Compute Scaling Parameters (Mean & SD)
mean_ELO_diff <- mean(df_long$ELO_diff, na.rm = TRUE)
sd_ELO_diff <- sd(df_long$ELO_diff, na.rm = TRUE)

mean_entry_diff <- mean(df_long$entry_diff, na.rm = TRUE)
sd_entry_diff <- sd(df_long$entry_diff, na.rm = TRUE)

# ✅ Create structured dataframe
scaling_results <- data.frame(
  Parameter = c("entry_diff_scaled_estimate", "mean_ELO_diff", "sd_ELO_diff", 
                "mean_entry_diff", "sd_entry_diff", 
                "max_entry_diff_threshold", "max_duration_threshold"),
  Value = c(entry_diff_scaled_estimate, mean_ELO_diff, sd_ELO_diff, 
            mean_entry_diff, sd_entry_diff, max_entry_diff, max_duration)
)

# ✅ Define Save Path with cohort_type in filename
save_filename <- paste0("scaling_parameters_", cohort_type, ".csv")
save_path <- file.path(data_dir, save_filename)

# ✅ Save as CSV
write_csv(scaling_results, save_path)

# ✅ Print Confirmation Message
message("Scaling parameters saved to: ", save_path)

# -------------------------------------------------------------
# 📈 9. Visualize Key Effects
# -------------------------------------------------------------

ggplot(df_long, aes(x = ELO_diff_scaled, y = Win)) +
  geom_smooth(method = "loess") +
  labs(title = "Winning Probability vs. ELO Difference",
       x = "ELO Difference (Standardized)",
       y = "Probability of Winning")

ggplot(df_long, aes(x = entry_diff_scaled, y = Win)) +
  geom_smooth(method = "loess") +
  labs(title = "Winning Probability vs. Entry Time Difference",
       x = "Entry Time Difference (Standardized)",
       y = "Probability of Winning")

ggplot(df_long, aes(x = duration_scaled, y = Win)) +
  geom_smooth(method = "loess") +
  labs(title = "Winning Probability vs. Own Duration",
       x = "Own Duration (Standardized)",
       y = "Probability of Winning")

# ✅ Histogram of entry_diff (after threshold-based NA replacement)
ggplot(df_long, aes(x = entry_diff)) +
  geom_histogram(binwidth = 0.1, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Histogram of entry_diff (With NA for Extreme Values)",
       x = "ELO Difference (Scaled)",
       y = "Frequency") +
  theme_minimal()

# -------------------------------------------------------------
# ✅ End of Script
# -------------------------------------------------------------

