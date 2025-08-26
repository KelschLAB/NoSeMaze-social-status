# Script to calculate characteristics of different NoSeMaze cohorts 
# (e.g., steepness, transitivity, correlation coefficients between ELO hierarchy ratings and David's score...)
# Jonathan Reinald, 17.06.2024

# Information:
# - based on aniDom and EloRating packages
# - input data is derived from the full_hierarchy.mat-files which we use in matlab and has to be transformed into sequence-format 
#   in a csv-file before using /zi-flstorage/group_entwbio/data/Shared/NoSeMaze/000_hierarchy/code for figures/functions/create_sequence_list.m
# - the script is currently adapted for cohort 1-10 from Sarah using day 1 to 14 (start and end day are adaptable)
# 
# Input: 
# - csv-files
# Output:
# 1. cohort_characteristics
# a data frame containing the cohort name, the optimal k, steepness, transitivity, linearity, stability and uncertainty for all cohorts
# 2. correlation_results
# a data frame containing the correlation results between (non-)randomized ELO ratings, David's score ...

# Load libraries (install these packages before install.packages("aniDom"))
library(aniDom)
library(EloRating)
library(dplyr)


# Pre-clearing
rm(list = ls())

# Predefine start and end day and k-value
start_day <- 1
end_day <- 14
K <- 100

# Specify the basic directory, of which the subfodlers contain the CSV files
base_dir <- "/myRootPath/NoSeMaze_Experiment_Social_Status/data/processed/" # --> CHANGE "myRootPath" MANUALLY TO YOUR ROOT DIRECTORY
results_dir <- "/myRootPath/NoSeMaze_Experiment_Social_Status/results/" # --> CHANGE "myRootPath" MANUALLY TO YOUR ROOT DIRECTORY

# Use list.files to recursively search for files matching the pattern
filelist <- list.files(
  path = base_dir,
  pattern = "^SequenceList_.*_TubeCompetitions\\.csv$", # Regex to match the file pattern
  full.names = TRUE, # Return full paths to files
  recursive = TRUE   # Search subdirectories
)

# Print the list of matching files
print(filelist)

# Optionally, save the list of files to a CSV for further inspection
write.csv(data.frame(FilePaths = filelist), 
          file = "filelist.csv", 
          row.names = FALSE)

# Initialize cohort_characteristics with correct lengths and types
cohort_characteristics <- data.frame(
  cohort = character(length(filelist)),         # Group names as character data
  optimal_K = numeric(length(filelist)),        # K as numeric data
  steepness = numeric(length(filelist)),        # Steepness as numeric data
  steepness_p = numeric(length(filelist)),      # Steepness significance as numeric data
  transitivity_tri = numeric(length(filelist)),     # Transitivity as numeric data
  transitivity_pt = numeric(length(filelist)),     # Transitivity as numeric data
  transitivity_p = numeric(length(filelist)),   # Transitivity significance as numeric data
  linearity_h1 = numeric(length(filelist)),        # Linearity as numeric data
  linearity_h2 = numeric(length(filelist)),        # Linearity as numeric data
  linearity_p = numeric(length(filelist)),      # Linearity significance as numeric data
  stabilityIndex = numeric(length(filelist)),   # StabilityIndex as numeric data
  numberofevents = numeric(length(filelist)),   # Number of events as numeric data
  uncertainty_split_mean = numeric(length(filelist)),
  uncertainty_split_CI_low = numeric(length(filelist)),
  uncertainty_split_CI_up = numeric(length(filelist)),
  uncertainty_rep = numeric(length(filelist)),
  stringsAsFactors = FALSE                      # Prevent automatic conversion to factors
)
# Add an empty list column for storing results
cohort_characteristics$res <- vector("list", length(filelist))
cohort_characteristics$InteractionMat <- vector("list", length(filelist))
cohort_characteristics$DavidsScore <- vector("list", length(filelist))
cohort_characteristics$DavidsScore_pij <- vector("list", length(filelist))
cohort_characteristics$DavidsScore_dij <- vector("list", length(filelist))
cohort_characteristics$sorted_non_randELOscore <- vector("list", length(filelist))
cohort_characteristics$sorted_randELOscore <- vector("list", length(filelist))
cohort_characteristics$sorted_DS <- vector("list", length(filelist))

# Initialize correlation_results frame with correct lengths and types
correlation_results <- data.frame(
  cohort = character(length(filelist)),         # Group names as character data
  randELOtoELO_estimate_p = numeric(length(filelist)),      # correlation coeff. between randELO and non-randELO as numeric data
  randELOtoELO_estimate_sp = numeric(length(filelist)),      # correlation coeff. between randELO and non-randELO as numeric data
  randELOtoELO_pval_p = numeric(length(filelist)),      # p-value between randELO and non-randELO as numeric data
  randELOtoELO_pval_sp = numeric(length(filelist)),      # p-value  between randELO and non-randELO as numeric data
  ELOtoDS_estimate_p = numeric(length(filelist)),      # correlation coeff. between DS and non-randELO as numeric data
  ELOtoDS_estimate_sp = numeric(length(filelist)),      # correlation coeff. between DS and non-randELO as numeric data
  ELOtoDS_pval_p = numeric(length(filelist)),      # p-value between DS and non-randELO as numeric data
  ELOtoDS_pval_sp = numeric(length(filelist)),      # p-value  between DS and non-randELO as numeric data
  randELOtoDS_estimate_p = numeric(length(filelist)),      # correlation coeff. between randELO and DS as numeric data
  randELOtoDS_estimate_sp = numeric(length(filelist)),      # correlation coeff. between randELO and DS as numeric data
  randELOtoDS_pval_p = numeric(length(filelist)),      # p-value between randELO and DS as numeric data
  randELOtoDS_pval_sp = numeric(length(filelist)),      # p-value  between randELO and DS as numeric data
  ELOaniDomtoELOEloRating_estimate_p = numeric(length(filelist)),      # correlation coeff. between randELO and DS as numeric data
  ELOaniDomtoELOEloRating_estimate_sp = numeric(length(filelist)),      # correlation coeff. between randELO and DS as numeric data
  ELOaniDomtoELOEloRating_pval_p = numeric(length(filelist)),      # p-value between randELO and DS as numeric data
  ELOaniDomtoELOEloRating_pval_sp = numeric(length(filelist)),      # p-value  between randELO and DS as numeric data  
  DSpijtoDSdij_estimate_p = numeric(length(filelist)),      # correlation coeff. between DSpij and DSdij as numeric data
  DSpijtoDSdij_estimate_sp = numeric(length(filelist)),      # correlation coeff. between DSpij and DSdij as numeric data
  DSpijtoDSdij_pval_p = numeric(length(filelist)),      # p-value between DSpij and DSdij as numeric data
  DSpijtoDSdij_pval_sp = numeric(length(filelist)),      # p-value  between DSpij and DSdij as numeric data 
  stringsAsFactors = FALSE                      # Prevent automatic conversion to factors
)

# Loop through each file and perform the Elo rating calculations as well as the calculation of the distinct characteristics
for (i in seq_along(filelist)) {
  file <- filelist[i]
  # Extract cohort name ("e.g., G1") from the file path
  cohort_name <- sub('.*SequenceList_(.*?)_days.*', '\\1', basename(file))
  
  # Initialize or update cohort_characteristics and correlation_results with the extracted cohort name
  cohort_characteristics$cohort[i] <- cohort_name
  correlation_results$cohort[i] <- cohort_name
  
  # Load the data
  data <- read.csv(file, header = TRUE)
  
  #  Extract winners and losers
  all_dates <- as.Date(data$day)
  unique_dates <- unique(all_dates)
  date_help <- which(all_dates == unique_dates[start_day])
  start_idx <- head(date_help,1)
  if (length(unique_dates) < end_day) {
    print("The length of unique dates is less than end_day. The last day of the recordings is used as last day.")
    # Perform your specific operation here
    date_help <- which(all_dates == tail(unique_dates,1))
    end_idx <- tail(date_help,1)
  } else {
    print("The length of unique dates is greater than or equal to end_day. The pre-selected end day is used as last day (as planned).")
    # Perform another operation or do nothing
    date_help <- which(all_dates == unique_dates[end_day])
    end_idx <- tail(date_help,1)
  }
  # restrict you data to the dates of interest
  winners <- data$winners[start_idx:end_idx]
  losers <- data$losers[start_idx:end_idx]
  date <- all_dates[start_idx:end_idx]
  cohort_characteristics$numberofevents[i] <- length(winners)
  
  ##############################################################################
  ################## EloRating package #########################################
  ##############################################################################
  
  # I. Calculate Elo ratings with the EloRating package
  
  # 1. k optimization
  res_prelimin <- elo.seq(winners, losers, date, draw = NULL, presence = NULL, startvalue = 0,
                          k = K, normprob = TRUE, init = "average", intensity = NULL,
                          iterate = 0, runcheck = TRUE, progressbar = FALSE)
  myCurrOptK <- optimizek(eloobject = res_prelimin, krange = c(2, 400), optimode = "loop", resolution = 200, itype = NULL, daterange = c(all_dates[start_idx],all_dates[end_idx]), 
                          burnin = 0, doplot = TRUE, progbar = TRUE)
  cohort_characteristics$optimal_K[i] <- myCurrOptK$best$k
  
  # 2. calculation of elo ratings
  cohort_characteristics$res[[i]] <- elo.seq(winners, losers, date, draw = NULL, presence = NULL, startvalue = 0,
                                             k = K, normprob = TRUE, init = "average", intensity = NULL,
                                             iterate = 0, runcheck = TRUE, progressbar = FALSE)
  # draw: logical, which interactions ended undecided
  # presence: optional data.frame, to supply data about presence and absence of individuals for part of the time the data collection covered. see details
  # startvalue: the value of Elo ratings of the two individuals that are involved in the first interaction of the overall sequence prior to this interaction. By default set to 1000.
  # k: factor k that determines the maximum change in ratings. By default k=100
  # normprob: logical (by default TRUE). Should a normal curve be assumed for calculating the winning/losing probablities, or a logistic curve. See winprob for details
  # init: character, what Elo rating does an individual have prior to its first interaction.Three options are available: average: individuals always start with the value specified in startvalue. Given stable composition of the group, this also reflects
  # the average Elo rating on each day in that group
  # intensity: intensity of interactions
  # iterate: not yet implemented
  
  # Extract elo ratings
  cat("Processing file", i, ":", file, "\n")
  elo_ratings <- extract_elo(cohort_characteristics$res[[i]], extractdate = all_dates[end_idx])
  print(elo_ratings)
  # extract_elo(eloobject = res, extractdate = date[15], standardize = FALSE, IDs = NULL, NA.interpolate = FALSE, daterange = 1)
  # -> NA.interpolate:
  # if FALSE (default), the last known rating is returned, which might not be from the
  # specified date itself (but older). If TRUE, ratings on days without observations are
  # linearly interpolated between days with known ratings (i.e. dates with observed interactions)
  # -> daterange:
  # if averaged ratings are desired, supply here the number of days from extractdate - 1. By default (daterange = 1), the ratings of the single extractdate
  # are returned. daterange = 2 produces average ratings from extractdate and the day after, and so on...
  
  # 2. plot (currently not working)
  # plt_filename <- paste0("eloplot_", i, ".png")
  # png(plot_filename)
  # eloplot(cohort_characteristics$res[[i]], ids="all", interpolate="yes", from="start", to="end", color=TRUE)
  
  # 3. Calculate interaction matrix (equivalent to the match_matrix in matlab)
  cohort_characteristics$InteractionMat[[i]]<- creatematrix(eloobject = cohort_characteristics$res[[i]], daterange = c(all_dates[start_idx],all_dates[end_idx]), drawmethod = "omit", onlyinteracting = FALSE, draw = NULL)
  
  # 4. Calculate DS
  cohort_characteristics$DavidsScore[[i]] <- DS(cohort_characteristics$InteractionMat[[i]], prop = c("Pij")) # Pij is what we use in matlab !!! "Dij" is the alternative option, 
  cohort_characteristics$DavidsScore_pij[[i]] <- DS(cohort_characteristics$InteractionMat[[i]], prop = c("Pij"))
  cohort_characteristics$DavidsScore_dij[[i]] <- DS(cohort_characteristics$InteractionMat[[i]], prop = c("Dij"))
  
  # 5. Steepness / transitivity / linearity / stability index
  mySteepness <- steepness(cohort_characteristics$InteractionMat[[i]], nrand = 10000, Dij = FALSE, returnfig = TRUE)
  myTransitivity <- transitivity(cohort_characteristics$InteractionMat[[i]], runs = 10000, returnfig = TRUE)
  myLinearity <- h.index(cohort_characteristics$InteractionMat[[i]], loops = 10000)
  myStabilityIndex <- stab_elo(cohort_characteristics$res[[i]], from = min(cohort_characteristics$res[[i]]$stability$date), to = all_dates[end_idx], weight = TRUE)
  cohort_characteristics$steepness[i] <- mySteepness[1]
  cohort_characteristics$steepness_p[i] <- mySteepness[3]
  cohort_characteristics$transitivity_pt[i] <- myTransitivity[1]
  cohort_characteristics$transitivity_tri[i] <- myTransitivity[2]
  cohort_characteristics$transitivity_p[i] <- myTransitivity[3]
  cohort_characteristics$linearity_h1[i] <- myLinearity$value[2]
  cohort_characteristics$linearity_h2[i] <- myLinearity$value[3]
  cohort_characteristics$linearity_p[i] <- myLinearity$value[5]
  cohort_characteristics$stabilityIndex[i] <- myStabilityIndex
  # Details
  # S ranges between 0 and 1, where 0 indicates an unstable hierarchy, in which the ordering reverses
  # every other day, and 1, in which the ordering is stable and no rank changes occur.
  # In contrast to the originally proposed S, this version is now standardized between 0 and 1, and
  # additionally, the interpretation is reversed, i.e. 1 refers to stable situations, whereas values closer to
  # 0 indicate more instable hierarchies
  
  ##############################################################################
  #################### aniDom package ##########################################
  ##############################################################################
  
  # 6. Calculation of randomised elo scores ("independent of order") with the package aniDom
  ELOscores <- elo_scores(winners, losers, identities = NULL, sigmoid.param = 1/100,
                          K = K, init.score = 0, randomise = TRUE, n.rands = 10000,
                          return.as.ranks = FALSE, return.trajectories = FALSE, dates = NULL)
  # !!! IMPORTANT: sigmoid.param --> to get similar results as in the EloRatingPackage, sigmoid.param has to be adapted (e.g. 1/200) !!!
  ## options
  # - sigmoid.param: A parameter of the Elo function that determines the steepness of the sigmoid function (i.e how much the scores change for small differences in rank). Smaller
  #   values flatten the shape (more linear), whereas larger values create a stronger threshold function (more change for small differences in rank)
  # - return.trajectories = TRUE and dates = date for trajectories
  # - randomise = TRUE and n.rands = 10000 if you want a randomised value, i.a., ignoring the order of the data (best comparable to David's score)
  # - optimal k estimation can be done with the EloRating package (see below)
  # randomised elo score is equivalent to the mean of the 10,000 randomisations
  mean_randomized_ELOscore <- rowMeans(ELOscores)
  
  # Convert the named numeric vector to a one-column matrix
  randomized_ELOscore <- matrix(mean_randomized_ELOscore, ncol = 1)
  # Set the row names
  rownames(randomized_ELOscore) <- names(mean_randomized_ELOscore)
  
  # save in cohort_characteristics
  cohort_characteristics$sorted_randELOscore[[i]] <- randomized_ELOscore
  
  # 7. Calculation of "normal" elo score
  non_randomized_ELOscore <- elo_scores(winners, losers, identities = NULL, sigmoid.param = 1/100,
                                        K = K, init.score = 0, randomise = FALSE, n.rands = 10000,
                                        return.as.ranks = FALSE, return.trajectories = FALSE, dates = NULL)
  non_randomized_ELOscore[,1] <- non_randomized_ELOscore[match(rownames(non_randomized_ELOscore), rownames(randomized_ELOscore))]
  rownames(non_randomized_ELOscore) <- rownames(non_randomized_ELOscore)[match(rownames(non_randomized_ELOscore), rownames(randomized_ELOscore))]
  
  # save in cohort_characteristics
  cohort_characteristics$sorted_non_randELOscore[[i]] <- non_randomized_ELOscore
  
  # sorted DS
  sortedDS <- matrix(cohort_characteristics$DavidsScore[[i]]$DS[match(rownames(randomized_ELOscore), cohort_characteristics$DavidsScore[[i]]$ID)], ncol = 1)
  rownames(sortedDS) <- cohort_characteristics$DavidsScore[[i]]$ID[match(rownames(randomized_ELOscore), cohort_characteristics$DavidsScore[[i]]$ID)]
  # save in cohort_characteristics
  cohort_characteristics$sorted_DS[[i]] <- sortedDS
  
  # sorted DS_dij
  sortedDS_dij <- matrix(cohort_characteristics$DavidsScore_dij[[i]]$DS[match(cohort_characteristics$DavidsScore_dij[[i]]$ID, cohort_characteristics$DavidsScore_pij[[i]]$ID)], ncol = 1)
  rownames(sortedDS_dij) <- cohort_characteristics$DavidsScore_dij[[i]]$ID[match(cohort_characteristics$DavidsScore_dij[[i]]$ID, cohort_characteristics$DavidsScore_pij[[i]]$ID)]
  
  
  # Check if all rownames are identical
  are_rownames_identical <- identical(rownames(sortedDS), rownames(non_randomized_ELOscore)) && 
    identical(rownames(non_randomized_ELOscore), rownames(randomized_ELOscore))
  # Print the result
  if (are_rownames_identical) {
    print("The rownames of sortedDS, non_randomized_ELOscore, and randomized_ELOscore are exactly the same.")
  } else {
    error_message <- "The rownames of sortedDS, non_randomized_ELOscore, and randomized_ELOscore are not the same."
  }
  
  # 8. Correlation between myDavidsScore and randomized_ELOscore
  # sorting and correlation between randomized ELO score and DS
  res_p <- cor.test(randomized_ELOscore,sortedDS,method = "pearson")
  res_sp <- cor.test(randomized_ELOscore,sortedDS,method = "spearman")
  correlation_results$randELOtoDS_estimate_p[i] <- res_p$estimate
  correlation_results$randELOtoDS_estimate_sp[i] <- res_sp$estimate
  correlation_results$randELOtoDS_pval_p[i] <- res_p$p.value
  correlation_results$randELOtoDS_pval_sp[i] <- res_sp$p.value  
  # correlation between non-randomized ELO score and DS
  res_p <- cor.test(non_randomized_ELOscore,sortedDS,method = "pearson")
  res_sp <- cor.test(non_randomized_ELOscore,sortedDS,method = "spearman")
  correlation_results$ELOtoDS_estimate_p[i] <- res_p$estimate
  correlation_results$ELOtoDS_estimate_sp[i] <- res_sp$estimate
  correlation_results$ELOtoDS_pval_p[i] <- res_p$p.value
  correlation_results$ELOtoDS_pval_sp[i] <- res_sp$p.value  
  # correlation between non-randomized ELO score and randomized ELO
  res_p <- cor.test(non_randomized_ELOscore,randomized_ELOscore,method = "pearson")
  res_sp <- cor.test(non_randomized_ELOscore,randomized_ELOscore,method = "spearman")
  correlation_results$randELOtoELO_estimate_p[i] <- res_p$estimate
  correlation_results$randELOtoELO_estimate_sp[i] <- res_sp$estimate
  correlation_results$randELOtoELO_pval_p[i] <- res_p$p.value
  correlation_results$randELOtoELO_pval_sp[i] <- res_sp$p.value    
  # correlation between DSdij and DSpij
  res_p <- cor.test(sortedDS_dij,cohort_characteristics$DavidsScore_pij[[i]]$DS,method = "pearson")
  res_sp <- cor.test(sortedDS_dij,cohort_characteristics$DavidsScore_pij[[i]]$DS,method = "spearman")
  correlation_results$DSpijtoDSdij_estimate_p[i] <- res_p$estimate
  correlation_results$DSpijtoDSdij_estimate_sp[i] <- res_sp$estimate
  correlation_results$DSpijtoDSdij_pval_p[i] <- res_p$p.value
  correlation_results$DSpijtoDSdij_pval_sp[i] <- res_sp$p.value    
  
  
  
  # correlation between non-randomized ELO score (aniDom) and ELO (EloRating)
  EloRating <- cohort_characteristics$res[[i]]$lmat[end_day-start_day+1,]
  sorted_EloRating <- EloRating[match(rownames(non_randomized_ELOscore),names(EloRating))]
  res_p <- cor.test(sorted_EloRating,non_randomized_ELOscore,method = "pearson")
  res_sp <- cor.test(sorted_EloRating,non_randomized_ELOscore,method = "spearman")
  correlation_results$ELOaniDomtoELOEloRating_estimate_p[i] <- res_p$estimate
  correlation_results$ELOaniDomtoELOEloRating_estimate_sp[i] <- res_sp$estimate
  correlation_results$ELOaniDomtoELOEloRating_pval_p[i] <- res_p$p.value
  correlation_results$ELOaniDomtoELOEloRating_pval_sp[i] <- res_sp$p.value  
  
  # 9. randomisation to estimate uncertainty (repeatability score)
  # Details
  # Each ordering of winners and losers will yield slightly different Elo scores. 
  # This function takes the Elo scores from n.rands randomisations of the order of interactions. 
  # It then computes the repeatability score. This repeatability score can provide some insight into the level of certainty (or robustness) of the input data. 
  # Our simulations suggest that a repeatability score above 0.8 suggests a reasonably robust hierarchy, given a large input dataset (can be unreliable for small datasets, i.e. < 10 observations per indiviudal), or for extremely flat hierarchies.
  cohort_characteristics$uncertainty_rep[i] <- estimate_uncertainty_by_repeatability(winners, losers, identities = NULL, 
                                                                                     sigmoid.param = 1/100, K = K, init.score = 0, n.rands = 10000)
  
  # 6. estimate uncertainty by splitting the dataset
  uncertainty_res <- estimate_uncertainty_by_splitting(winners, losers, identities = NULL,
                                                       sigmoid.param = 1/100, K = K, init.score = 0, randomise = TRUE, n.rands = 10000)
  cohort_characteristics$uncertainty_split_mean[i] <- uncertainty_res[1]
  cohort_characteristics$uncertainty_split_CI_low[i] <- uncertainty_res[2]
  cohort_characteristics$uncertainty_split_CI_up[i] <- uncertainty_res[3]
}


# Save as .RData files as in your original code
save(correlation_results, file = file.path(results_dir, paste0("CorrelationResults_", start_day, "_to_", end_day, ".RData")))
save(cohort_characteristics, file = file.path(results_dir, paste0("CohortCharacteristics_", start_day, "_to_", end_day, ".RData")))

# Flattening list columns
cohort_characteristics <- cohort_characteristics %>%
  mutate(across(where(is.list), ~sapply(., function(x) if (length(x) == 1) x else paste(x, collapse = ";"))))
# Now save as CSV
write.csv(cohort_characteristics, file = file.path(results_dir, paste0("CohortCharacteristics_", start_day, "_to_", end_day, ".csv")), row.names = FALSE)

# Flattening list columns
correlation_results <- correlation_results %>%
  mutate(across(where(is.list), ~sapply(., function(x) if (length(x) == 1) x else paste(x, collapse = ";"))))
# Now save as CSV
write.csv(correlation_results, file = file.path(results_dir, paste0("CorrelationResults_", start_day, "_to_", end_day, ".csv")), row.names = FALSE)

