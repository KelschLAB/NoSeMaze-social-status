# Source script for cohort calculations
# Author: Jonathan Reinald, Revised: 17.06.2024

# Libraries
library(aniDom)
library(EloRating)
library(dplyr)

# Function to filter files based on date range
filter_files_by_date_range <- function(filelist, day_range_name, cohorts_tbl) {
  filtered_files <- vector("character")
  for (file in filelist) {
    data <- read.csv(file, header = TRUE)
    
    # Extract cohort identifier from the file name (assuming it contains cohort info)
    cohort_name <- gsub(".*_(cohort\\d+)_.*", "\\1", basename(file)) # Adjust regex if needed
    
    # Find the respective start_date for the cohort
    cohort_row <- cohorts_tbl[cohorts_tbl$cohort == cohort_name, ]
    if (nrow(cohort_row) == 0) {
      warning(sprintf("Cohort not found for file: %s. Skipping.", file))
      next
    }
    cohort_start_date <- as.Date(cohort_row$start_date, format = "%d.%m.%Y")
    num_days <- as.numeric(cohort_row$num_days)
    day_range <- extract_day_range(day_range_name, num_days)
    if (is.null(day_range)) next
    
    start_range <- day_range[1]
    end_range <- day_range[2]
    
    # Filter dates within the desired range
    valid_dates <- unique(as.Date(data$day))
    valid_dates <- valid_dates[valid_dates >= cohort_start_date + (start_range - 1) & 
                                 valid_dates <= cohort_start_date + (end_range - 1)]
    
    # Calculate total days in the range
    total_days <- as.numeric(cohort_start_date + (end_range - 1) - (cohort_start_date + (start_range - 1)) + 1)
    
    # Check if data is available for more than 50% of the days in the range
    if (length(valid_dates) >= 1) {
      filtered_files <- c(filtered_files, file)
    }
  }
  
  return(filtered_files)
}

# Function to extract numeric day ranges from MATLAB-style labels
extract_day_range <- function(range_name, num_days) {
  if (grepl("End", range_name)) {
    start_day <- as.numeric(sub("D(\\d+)_End", "\\1", range_name))
    end_day <- num_days
  } else if (range_name == "Last14") {
    start_day <- max(num_days - 13, 1)  # Ensure we don't go before day 1
    end_day <- num_days
  } else if (grepl("D\\d+_\\d+", range_name)) {
    matches <- regmatches(range_name, gregexpr("\\d+", range_name))[[1]]
    start_day <- as.numeric(matches[1])
    end_day <- min(as.numeric(matches[2]), num_days)
  } else {
    return(NULL)  # Skip invalid ranges
  }
  return(c(start_day, end_day))
}


# Function to calculate cohort characteristics
calculate_cohort_characteristics <- function(filelist, day_range_name, K, cohorts_tbl) {
  # Initialize data frames and list columns
  cohort_characteristics <- data.frame(
    cohort = character(length(filelist)),
    numberofevents = numeric(length(filelist)),
    optimal_K = numeric(length(filelist)),
    steepness = numeric(length(filelist)),
    steepness_p = numeric(length(filelist)),
    transitivity_tri = numeric(length(filelist)),
    transitivity_pt = numeric(length(filelist)),
    transitivity_p = numeric(length(filelist)),
    linearity_h1 = numeric(length(filelist)),
    linearity_h2 = numeric(length(filelist)),
    linearity_p = numeric(length(filelist)),
    stabilityIndex = numeric(length(filelist)),
    uncertainty_split_mean = numeric(length(filelist)),
    uncertainty_split_CI_low = numeric(length(filelist)),
    uncertainty_split_CI_up = numeric(length(filelist)),
    uncertainty_rep = numeric(length(filelist)),
    stringsAsFactors = FALSE
  )
  
  # Add list columns for complex objects
  cohort_characteristics$res <- vector("list", length(filelist))
  cohort_characteristics$InteractionMat <- vector("list", length(filelist))
  cohort_characteristics$DavidsScore <- vector("list", length(filelist))
  cohort_characteristics$sorted_randELOscore <- vector("list", length(filelist))
  cohort_characteristics$sorted_non_randELOscore <- vector("list", length(filelist))
  cohort_characteristics$sorted_DS <- vector("list", length(filelist))
  cohort_characteristics$sorted_IDs <- vector("list", length(filelist))
  
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
  
  # Iterate over each file
  for (i in seq_along(filelist)) {

    # Current file
    file <- filelist[i]
    
    # Extract cohort name
    cohort_name <- sub('.*SequenceList_(.*?)_days.*', '\\1', basename(file))
    cohort_characteristics$cohort[i] <- cohort_name
    correlation_results$cohort[i] <- cohort_name
    
    # Find the respective start_date for the cohort
    cohort_row <- cohorts_tbl[cohorts_tbl$cohort == cohort_name, ]
    if (nrow(cohort_row) == 0) {
      warning(sprintf("Cohort not found for file: %s. Skipping.", file))
      next
    }
    cohort_start_date <- as.Date(cohort_row$start_date, format = "%d.%m.%Y")
    num_days <- as.numeric(cohort_row$num_days)
    day_range <- extract_day_range(day_range_name, num_days)
    if (is.null(day_range)) next
    
    start_range <- day_range[1]
    end_range <- day_range[2]
    
    # Load data
    data <- read.csv(file, header = TRUE)
    all_dates <- as.Date(data$day)
    
    # Number of days in the range
    n_day_in_range <- end_range - start_range + 1
    
    # Determine start and end indices
    start_idx <- which(all_dates == as.Date(cohort_start_date + (start_range - 1)))[1]
    end_idx <- tail(which(all_dates == as.Date(cohort_start_date + (end_range - 1))), 1)
    
    # Iteratively adjust start_idx if empty
    if (length(start_idx) == 0) {
      for (offset in 1:n_day_in_range) {
        start_idx <- which(all_dates == as.Date(cohort_start_date + (start_range - 1 + offset)))[1]
        if (length(start_idx) > 0) break
      }
    }
    
    # Iteratively adjust end_idx if empty
    if (length(end_idx) == 0) {
      for (offset in 1:n_day_in_range) {
        end_idx <- tail(which(all_dates == as.Date(cohort_start_date + (end_range - 1 - offset))), 1)
        if (length(end_idx) > 0) break
      }
    }

    # Final check: Skip if indices are still empty
    if (length(start_idx) == 0 || length(end_idx) == 0) {
      warning(sprintf("Could not determine valid indices for cohort starting on %s in the given range.", cohort_start_date))
      next
    }
    
    # Restrict data to the dates of interest
    winners <- data$winners[start_idx:end_idx]
    losers <- data$losers[start_idx:end_idx]
    date <- all_dates[start_idx:end_idx]
    cohort_characteristics$numberofevents[i] <- length(winners)
    
    # Check for sufficient events per day
    total_days <- length(unique(date))
    total_events <- length(winners)
    if (total_days > 0 && (total_events / total_days) < 3) {
      warning(paste("Cohort", cohort_name, "has fewer than 3 events per day on average."))
    }
    
    ##############################################################################
    ################## EloRating package #########################################
    ##############################################################################
    
    # I. Calculate Elo ratings with the EloRating package
    # 1. k optimization
    res_prelimin <- elo.seq(winners, losers, date, draw = NULL, presence = NULL, startvalue = 0,
                            k = 100, normprob = TRUE, init = "average", intensity = NULL,
                            iterate = 0, runcheck = TRUE, progressbar = FALSE)
    myCurrOptK <- optimizek(eloobject = res_prelimin, krange = c(2, 400), optimode = "loop", resolution = 200, itype = NULL, daterange = c(all_dates[start_idx],all_dates[end_idx]), 
                            burnin = 0, doplot = TRUE, progbar = TRUE)
    cohort_characteristics$optimal_K[i] <- myCurrOptK$best$k
    
    # 2. calculation of elo ratings
    if (is.character(K) && K == 'optimal') {
      cohort_characteristics$res[[i]] <- elo.seq(winners, losers, date, draw = NULL, presence = NULL, startvalue = 0,
                                                 k = cohort_characteristics$optimal_K[i], normprob = TRUE, init = "average", intensity = NULL,
                                                 iterate = 0, runcheck = TRUE, progressbar = FALSE)
    } else {
      cohort_characteristics$res[[i]] <- elo.seq(winners, losers, date, draw = NULL, presence = NULL, startvalue = 0,
                                                 k = K, normprob = TRUE, init = "average", intensity = NULL,
                                                 iterate = 0, runcheck = TRUE, progressbar = FALSE)
    }
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
    
    # 3. Calculate interaction matrix (equivalent to the match_matrix in matlab)
    cohort_characteristics$InteractionMat[[i]]<- creatematrix(eloobject = cohort_characteristics$res[[i]], daterange = c(all_dates[start_idx],all_dates[end_idx]), drawmethod = "omit", onlyinteracting = FALSE, draw = NULL)
    
    # 4. Calculate DS
    cohort_characteristics$DavidsScore[[i]] <- DS(cohort_characteristics$InteractionMat[[i]], prop = c("Pij")) # Pij is what we use in matlab !!! "Dij" is the alternative option, 
    cohort_characteristics$DavidsScore_pij[[i]] <- DS(cohort_characteristics$InteractionMat[[i]], prop = c("Pij"))
    cohort_characteristics$DavidsScore_dij[[i]] <- DS(cohort_characteristics$InteractionMat[[i]], prop = c("Dij"))
    
    # 5. Steepness / transitivity / linearity / stability index
    steepness_res <- steepness(cohort_characteristics$InteractionMat[[i]], nrand = 10000, Dij = FALSE, returnfig = FALSE)
    cohort_characteristics$steepness[i] <- steepness_res[1]
    cohort_characteristics$steepness_p[i] <- steepness_res[3]
    
    transitivity_res <- transitivity(cohort_characteristics$InteractionMat[[i]], runs = 10000, returnfig = FALSE)
    cohort_characteristics$transitivity_pt[i] <- transitivity_res[1]
    cohort_characteristics$transitivity_tri[i] <- transitivity_res[2]
    cohort_characteristics$transitivity_p[i] <- transitivity_res[3]
    
    linearity_res <- h.index(cohort_characteristics$InteractionMat[[i]], loops = 10000)
    cohort_characteristics$linearity_h1[i] <- linearity_res$value[2]
    cohort_characteristics$linearity_h2[i] <- linearity_res$value[3]
    cohort_characteristics$linearity_p[i] <- linearity_res$value[5]
    
    stability_res <- stab_elo(cohort_characteristics$res[[i]], from = min(cohort_characteristics$res[[i]]$stability$date), to = all_dates[end_idx], weight = TRUE)
    cohort_characteristics$stabilityIndex[i] <- stability_res
    
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
    if (is.character(K) && K == 'optimal') {
      ELOscores <- elo_scores(winners, losers, identities = NULL, sigmoid.param = 1/100,
                              K = cohort_characteristics$optimal_K[i], init.score = 0, randomise = TRUE, n.rands = 1000,
                              return.as.ranks = FALSE, return.trajectories = FALSE, dates = NULL)
    } else {
      ELOscores <- elo_scores(winners, losers, identities = NULL, sigmoid.param = 1/100,
                              K = K, init.score = 0, randomise = TRUE, n.rands = 1000,
                              return.as.ranks = FALSE, return.trajectories = FALSE, dates = NULL)
    }

    # !!! IMPORTANT: sigmoid.param --> to get similar results as in the EloRatingPackage, sigmoid.param has to be adapted (e.g. 1/200) !!!
    ## options
    # - sigmoid.param: A parameter of the Elo function that determines the steepness of the sigmoid function (i.e how much the scores change for small differences in rank). Smaller
    #   values flatten the shape (more linear), whereas larger values create a stronger threshold function (more change for small differences in rank)
    # - return.trajectories = TRUE and dates = date for trajectories
    # - randomise = TRUE and n.rands = 1000 if you want a randomised value, i.a., ignoring the order of the data (best comparable to David's score)
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
    if (is.character(K) && K == 'optimal') {
      non_randomized_ELOscore <- elo_scores(winners, losers, identities = NULL, sigmoid.param = 1/100,
                                            K = cohort_characteristics$optimal_K[i], init.score = 0, randomise = FALSE, n.rands = 1000,
                                            return.as.ranks = FALSE, return.trajectories = FALSE, dates = NULL)
    } else {
      non_randomized_ELOscore <- elo_scores(winners, losers, identities = NULL, sigmoid.param = 1/100,
                                            K = K, init.score = 0, randomise = FALSE, n.rands = 1000,
                                            return.as.ranks = FALSE, return.trajectories = FALSE, dates = NULL)
    }
    non_randomized_ELOscore[,1] <- non_randomized_ELOscore[match(rownames(non_randomized_ELOscore), rownames(randomized_ELOscore))]
    rownames(non_randomized_ELOscore) <- rownames(non_randomized_ELOscore)[match(rownames(non_randomized_ELOscore), rownames(randomized_ELOscore))]
    
    # save in cohort_characteristics
    cohort_characteristics$sorted_non_randELOscore[[i]] <- non_randomized_ELOscore
    
    # sorted DS
    sortedDS <- matrix(cohort_characteristics$DavidsScore[[i]]$DS[match(rownames(randomized_ELOscore), cohort_characteristics$DavidsScore[[i]]$ID)], ncol = 1)
    rownames(sortedDS) <- cohort_characteristics$DavidsScore[[i]]$ID[match(rownames(randomized_ELOscore), cohort_characteristics$DavidsScore[[i]]$ID)]
    # save in cohort_characteristics
    cohort_characteristics$sorted_DS[[i]] <- sortedDS
    # add sorted IDs
    cohort_characteristics$sorted_IDs[[i]] <- rownames(sortedDS)
    
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
    # Access the final row of the lmat matrix
    EloRating <- cohort_characteristics$res[[i]]$lmat[nrow(cohort_characteristics$res[[i]]$lmat), ]
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
    if (is.character(K) && K == 'optimal') {
      cohort_characteristics$uncertainty_rep[i] <- estimate_uncertainty_by_repeatability(winners, losers, identities = NULL, 
                                                                                         sigmoid.param = 1/100, K = cohort_characteristics$optimal_K[i], init.score = 0, n.rands = 1000)
    } else {
      cohort_characteristics$uncertainty_rep[i] <- estimate_uncertainty_by_repeatability(winners, losers, identities = NULL, 
                                                                                         sigmoid.param = 1/100, K = K, init.score = 0, n.rands = 1000)
    }
    
    # 6. estimate uncertainty by splitting the dataset
    if (is.character(K) && K == 'optimal') {
      uncertainty_res <- estimate_uncertainty_by_splitting(winners, losers, identities = NULL,
                                                           sigmoid.param = 1/100, K = cohort_characteristics$optimal_K[i], init.score = 0, randomise = TRUE, n.rands = 1000)
    } else {
      uncertainty_res <- estimate_uncertainty_by_splitting(winners, losers, identities = NULL,
                                                           sigmoid.param = 1/100, K = K, init.score = 0, randomise = TRUE, n.rands = 1000)
    }

    cohort_characteristics$uncertainty_split_mean[i] <- uncertainty_res[1]
    cohort_characteristics$uncertainty_split_CI_low[i] <- uncertainty_res[2]
    cohort_characteristics$uncertainty_split_CI_up[i] <- uncertainty_res[3]
    browser()
  }
  print(list(cohort_characteristics = cohort_characteristics, correlation_results = correlation_results))
  return(list(cohort_characteristics = cohort_characteristics, correlation_results = correlation_results))
}



