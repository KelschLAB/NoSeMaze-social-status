README – Adjusting ELO Ratings for entry time difference in tube competitions
-------------------------------------------------------------------------------------------------------
This directory contains scripts for merging, LME-modeling, and adjusting ELO ratings from tube competition data in the NoSeMaze experiment.
SOPs are saved in: myRootPath\NoSeMaze_Experiment_Socia_Status\docs\SOPs\adjusted_ELO
-------------------------------------------------------------------------------------------------------
A. 
Script Overview & Purpose


📌 Step 1: Merge ELO Ratings with Competition Data
Script:
	📂 analysis/scripts/adjusted_ELO/combine_sequence_files_with_ELO_ratings.R
Purpose:
- Merges stepwise ELO ratings (pre/post-competition) with tube competition sequence files (entry/exit times, durations).
- Prepares structured data for linear mixed modeling (LME).
- Either use K100 for all cohorts or Koptimal (optimal K-value of each cohort) for the ELO ratings
- Saves merged dataset as:
	📂 data/processed/cross_cohort_files/tube/merged_elo_sequence_data/Merged_Tube_Competitions_Koptimal.csv


📌 Step 2: Analyze Tube Competition Outcomes Using LME
Script:
	📂 analysis/scripts/adjusted_ELO/LME_model_tube_competition_outcome.R
Purpose:
- Uses Generalized Linear Mixed Models (GLMMs) to assess the influence of ELO difference and entry time on winning probability.
- Computes scaling parameters for ELO adjustments.
- Saves scaling parameters as:
	📂 data/processed/cross_cohort_files/tube/merged_elo_sequence_data/scaling_parameters_Koptimal.csv


📌 Step 3: Calculate Adjusted ELO Ratings
Script:
	📂 analysis/scripts/adjusted_ELO/calculate_adjusted_ELO_ratings.R
Purpose:
- Computes standard and corrected ELO ratings for tube competitions. Script needs the scaling parameters from the step before.
- IMPORTANT: Adjusts ELO updates based on entry time differences in a step-wise manner.
- Saves results as:
	📂 data/processed/cross_cohort_files/tube/merged_elo_sequence_data/Final_ELO_ratings_by_cohort_Koptimal.csv
- Generates & saves correlation plots between corrected & uncorrected ELO:
	📂 results/figures/cross_cohort/adjusted_for_entry_difference/ELO_Correlation_Corrected_vs_Uncorrected_Koptimal.pdf




B.
Data Paths & Output Files

📂 Input Data:
- data/processed/cross_cohort_files/tube/
	- CohortCharacteristics_TubeCompetitions_D1_End_Koptimal.Rdata (Required for ELO calculations)
	- Merged_Tube_Competitions_Koptimal.csv (Merged dataset for LME analysis)
📂 Output Data:
- data/processed/cross_cohort_files/tube/merged_elo_sequence_data/
	- Final_ELO_ratings_by_cohort_Koptimal.csv (Final computed ELO ratings)
	- scaling_parameters_Koptimal.csv (Scaling values for modeling)
📂 Figures & Results:
- results/figures/cross_cohort/adjusted_for_entry_difference/
	- ELO_Correlation_Corrected_vs_Uncorrected_Koptimal.pdf (Correlation plot between standard & corrected ELO)


C.
Execution Order

1.Run: combine_sequence_files_with_ELO_ratings.R → Merges ELO ratings with competition data containing entry information.
2️.Run: LME_model_competition_outcome.R → Performs LME analysis & saves scaling parameters for adjustment of ELO ratings in the next step.
3️.Run: calculate_adjusted_ELO_ratings.R → Computes standard & corrected ELO.

