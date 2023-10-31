




#' process_a_study
#'
#' Main function for calculating time series features for study parameters, and flagging study sites with systematic bias.
#'
#' @param subjects Data frame with one row per study subject.
#' @param parameters Data frame with one row per study parameter (e.g. labs and vital signs) to process.
#' @param data Measurements of study parameters. One row per data point.
#' @param custom_timeseries Pre-defined time series.
#' @param timeseries_features_to_calculate Vector of time series features to calculate.
#' @param default_minimum_timepoints_per_series Default minimum timepoints per time series. This is used if no parameter-specific minimum is defined in parameters.
#' @param default_minimum_subjects_per_series Default minimum of subjects per time series. This is used if no parameter-specific minimum is defined in parameters.
#' @param default_max_share_missing_timepoints_per_series Maximum share of time points which can be missing per subject. This is used if no parameter-specific minimum is defined in parameters.
#' @param default_generate_change_from_baseline If set TRUE, time series and their features are calculate also for change-from-baseline values.
#' @param autogenerate_timeseries If set to TRUE, automatic definition of time series is used. If set to FALSE, custom_timeseries must have at least one time series defined.
#' @return List with four data frames. Timeseries: definition of the time series used. Timeseries_features: features calculated from the time series. PCA_coordinates: principal components of individual time series for visualizing similarity. Site_scores: biasness scores for sites.
#'
#' @export
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
process_a_study <- function(subjects, parameters, data, custom_timeseries, timeseries_features_to_calculate,
                            default_minimum_timepoints_per_series, default_minimum_subjects_per_series,
                            default_max_share_missing_timepoints_per_series, default_generate_change_from_baseline,
                            autogenerate_timeseries) {

  # Disable confusing "`summarise()` has grouped output by 'XX'. You can override using the `.groups` argument." messages.
  options(dplyr.summarise.inform = FALSE)

  # Makes sure all the necessary input is there and that the data is in correct format etc.
  check_input_data(subjects, parameters, data, custom_timeseries, timeseries_features_to_calculate,
                   default_minimum_timepoints_per_series, default_minimum_subjects_per_series,
                   default_max_share_missing_timepoints_per_series, default_generate_change_from_baseline,
                   autogenerate_timeseries)

  # Replace empty parameter settings with defaults.
  parameters <- parameters %>%
    mutate(time_point_count_min = ifelse(is.na(.data$time_point_count_min), .env$default_minimum_timepoints_per_series, .data$time_point_count_min),
           subject_count_min = ifelse(is.na(.data$subject_count_min), .env$default_minimum_subjects_per_series, .data$subject_count_min),
           max_share_missing = ifelse(is.na(.data$max_share_missing), .env$default_max_share_missing_timepoints_per_series, .data$max_share_missing),
           generate_change_from_baseline = ifelse(is.na(.data$generate_change_from_baseline), .env$default_generate_change_from_baseline, .data$generate_change_from_baseline)
    )

  # This mapping table is used later for mapping parameter-specific time ranks into human-friendly names.
  timepoint_rank_to_name_mapping <- data %>%
    mutate(names = paste0(.data$timepoint_1_name, "_", .data$timepoint_2_name)) %>%
    distinct(.data$parameter_id, .data$timepoint_rank, names)

  # Generate a table of subjects and the time points they have data for.
  timepoints_and_subjects <- data %>%
    filter(!is.na(.data$result) & .data$result != '') %>%
    mutate(has_baseline_value = ifelse(!is.na(.data$baseline) & .data$baseline != '', "Yes", "No")) %>%
    arrange(.data$parameter_id, .data$subject_id, .data$timepoint_rank, .data$has_baseline_value) %>%
    distinct(.data$parameter_id, .data$subject_id, .data$timepoint_rank, .data$has_baseline_value)

  if(autogenerate_timeseries) {

    # Autogenerate the time point combos (time series) for each parameter.
    timeseries_per_param <- timepoints_and_subjects %>%
      group_by(.data$parameter_id) %>%
      nest() %>%
      left_join(y=select(parameters, c("parameter_id", "time_point_count_min", "subject_count_min", "max_share_missing", "generate_change_from_baseline")), by="parameter_id") %>%
      rename(dataset = .data$data) %>%
      mutate(baseline = ifelse(.data$generate_change_from_baseline == TRUE, list(c("original", "cfb")), "original")) %>%
      unnest("baseline") %>%
      rowwise() %>%
      mutate(output = list(pick_timepoint_combos(.data$dataset, .data$time_point_count_min, .data$subject_count_min, .data$max_share_missing, .data$baseline))) %>%
      unnest("output") %>%
      ungroup() %>%
      mutate(timeseries_id = paste0('ts_', row_number(), '_autogen')) %>%
      mutate(timeseries_id = paste0(.data$timeseries_id, '_', .data$baseline)) %>%
      select(- c("dataset", "time_point_count_min", "subject_count_min", "generate_change_from_baseline", "max_share_missing"))

    # If one of the custom timeseries was generated automatically, remove it from the autogenerated time series list
    # and add it later together with all the custom time series.
    if( nrow(custom_timeseries) > 0 ) {

      timeseries_per_param <- timeseries_per_param %>%
        anti_join(y=select(custom_timeseries, c("parameter_id", "timepoint_combo")), by = NULL) # Remove any autogenerated timeseries that are in the list of custom timeseries

    }

  }

  # Process any custom time series. Mainly this checks whether enough subjects have enough data for the time series.
  if( nrow(custom_timeseries) > 0 ) {

    timeseries_per_param_custom <- custom_timeseries %>%
      left_join(y=select(parameters, c("parameter_id", "subject_count_min", "max_share_missing", "generate_change_from_baseline")), by="parameter_id") %>%
      rowwise() %>%
      mutate(baseline = ifelse(.data$generate_change_from_baseline == TRUE, list(c("original", "cfb")), "original")) %>%
      unnest("baseline") %>%
      mutate(timepoint_combo_subjects = pick_subjects_for_custom_timeseries(.data$timepoints_and_subjects, .data$timepoint_combo, .data$max_share_missing, .data$parameter_id, .data$baseline)) %>%
      filter(str_count(.data$timepoint_combo_subjects, ";") + 1 >= .data$subject_count_min) %>% # Only pass time series with enough subjects.
      mutate(timeseries_id = paste0(.data$timeseries_id, '_', .data$baseline)) %>%
      select(-c("subject_count_min", "generate_change_from_baseline", "max_share_missing"))

    if(autogenerate_timeseries) {

      # Merge timeseries tables.
      timeseries_per_param <- bind_rows(timeseries_per_param, timeseries_per_param_custom)

    } else {

      timeseries_per_param <- timeseries_per_param_custom

    }

  }



  # If there are no time series to process, stop
  if(nrow(timeseries_per_param) == 0) {

    return(list("timeseries" = NULL, "timeseries_features" = NULL,
                "PCA_coordinates" = NULL, "site_scores" = NULL))

  }

  # Count the number of time points in each time series
  timeseries_per_param$timepoint_count <- str_count(timeseries_per_param$timepoint_combo, ";") + 1


  # Calculate time series features and the principal components used in similarity plots.
  results <- timeseries_per_param %>%
    filter(.data$baseline == "original" | .data$timepoint_count > 1) %>% # The latter condition requires that change-from-baseline series must have more than one time point.
    rowwise() %>%
    mutate(timeseries_wide = list(generate_wide_timeseries_table(.data$parameter_id, .data$timepoint_combo, .data$timepoint_combo_subjects, .data$baseline, .env$data))) %>%
    filter(nrow(.data$timeseries_wide) > 0) %>%
    mutate(principal_components = list(calculate_principal_components(.data$timeseries_wide))) %>%
    mutate(ts_features = list(calculate_ts_features(.data$timeseries_wide, .data$baseline, .env$timeseries_features_to_calculate, .env$subjects))) %>%
    ungroup()

  # Parse the timepoint combination string into a more readable format.
  readable_timepoint_combo_names <- results %>%
    distinct(.data$timepoint_combo, .data$parameter_id) %>%
    rowwise() %>%
    mutate(timepoint_combo_readable = parse_readable_timeseries_combo_string(.data$timepoint_combo, .data$parameter_id, .env$timepoint_rank_to_name_mapping))

  # Generate a list of time series evaluated.
  tso_timeseries <- results %>%
    left_join(y=readable_timepoint_combo_names, by=c("parameter_id", "timepoint_combo")) %>%
    select(c("timeseries_id", "parameter_id", "baseline", "timepoint_combo", "timepoint_combo_readable", "timepoint_count"))

  # Parse time series features into a proper format.
  tso_features <- results %>%
    select(c("timeseries_id", "ts_features")) %>%
    unnest("ts_features") %>%
    left_join(y=subjects, by="subject_id") %>%
    rename(c(feature_value = "value")) %>%
    relocate(c("timeseries_id", "subject_id", "feature", "feature_value"))

  # Parse the PCA coordinates into its own result table.
  tso_pca_coordinates <- results %>%
    select(c("timeseries_id", "principal_components")) %>%
    unnest("principal_components") %>%
    left_join(y=subjects, by="subject_id") %>%
    rename(c(pc1 = "PC1", pc2 = "PC2")) %>%
    select(c("timeseries_id", "subject_id", "pc1", "pc2"))

  # Calculate the biasness scores for sites
  tso_site_scores <- results %>%
    select(-c("timeseries_wide", "principal_components")) %>%
    unnest("ts_features") %>%
    inner_join(y=subjects, by="subject_id") %>%
    left_join(y=select(parameters, c("parameter_id", "subject_count_min")), by="parameter_id") %>% # Look up minimum subjects expected per parameter.
    group_by(.data$timeseries_id, .data$feature) %>%
    filter( n_distinct(.data$site) >= 2 & n() >= .data$subject_count_min) %>%
    nest()


  if(nrow(tso_site_scores) >= 1) {

    # Process the site scores further only if there is something to process...

    tso_site_scores <- tso_site_scores %>%
      rowwise() %>%
      mutate(site_p_values = list(calculate_site_bias_ts_features(.data$feature, .data$data))) %>%
      select(-data) %>%
      unnest("site_p_values") %>%
      mutate(pvalue_kstest = as.numeric(.data$pvalue_kstest), kstest_statistic = as.numeric(.data$kstest_statistic)) %>%
      ungroup() %>%
      mutate(fdr_adjusted_pvalue_ks = p.adjust(.data$pvalue_kstest, method = "fdr"),
             reference_group = "all") %>% # Correct for multiple testing
      mutate(pvalue_kstest_logp = -log10(.data$pvalue_kstest), fdr_corrected_pvalue_logp = -log10(.data$fdr_adjusted_pvalue_ks)) %>%
      mutate(pvalue_kstest_logp = if_else( is.infinite(.data$pvalue_kstest_logp), 30, .data$pvalue_kstest_logp),
             fdr_corrected_pvalue_logp = if_else( is.infinite(.data$fdr_corrected_pvalue_logp), 30, .data$fdr_corrected_pvalue_logp)) %>%
      select(c("timeseries_id", "site", "country", "feature", "pvalue_kstest_logp", "kstest_statistic", "fdr_corrected_pvalue_logp", "reference_group", "subj_count")) %>%
      rename(c(subject_count = "subj_count"))

  } else {

    # Otherwise create an empty data frame.

    tso_site_scores <- data.frame(matrix(ncol = 9, nrow = 0))
    colnames(tso_site_scores) <- c("timeseries_id", "site", "country", "feature", "pvalue_kstest_logp", "kstest_statistic", "fdr_corrected_pvalue_logp", "reference_group", "subject_count")

  }

  # Combine the four result tables into a list which the function returns.
  return(list("timeseries" = tso_timeseries, "timeseries_features" = tso_features,
              "PCA_coordinates" = tso_pca_coordinates, "site_scores" = tso_site_scores))

}



#' calculate_site_bias_ts_features
#'
#' Performs Kolmogorov-Smirnov (KS) test for a particular parameter and time series feature.
#'
#' @param this_feature Type of the time series feature.
#' @param this_data Feature values for each subject.
#' @return Data frame with KS p-values for each site.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
calculate_site_bias_ts_features <- function(this_feature, this_data) {

  # Add a tiny random number to each value to avoid problems ks.test has with ties.
  this_data$value <- this_data$value + rnorm(nrow(this_data), mean = 0, sd = 0.00001)

  site_pvalues <- data.frame()

  # Define the hypothesis (i.e. whether we are looking at one-sided or two-sided bias).
  nullhypo_ks <- case_when(
    this_feature == "own_site_simil_score" ~ "less",
    this_feature == "outlier_score" ~ "less",
    this_feature == "unique_value_count_relative" ~ "greater",
    TRUE ~ "two.sided"
  )

  site_metadata <- this_data %>%
    group_by(.data$site) %>%
    summarise(country = first(.data$country), subj_count = n_distinct(.data$subject_id))

  # Put site name and value columns into vectors. This makes the for loop much faster.
  unique_sites <- unique(this_data$site)
  sites <- this_data$site
  values <- this_data$value

  for(this_site in unique_sites) {

    this_site_value_indices <- which(sites == this_site)

    within_set_values <- values[this_site_value_indices]
    outside_set_values <- values[-this_site_value_indices]

    ks_results <- ks.test(x = within_set_values, y = outside_set_values, alternative = nullhypo_ks)

    site_pvalues <- bind_rows(site_pvalues,
                              c("site" = this_site,
                                "pvalue_kstest" = ks_results$p.value, "kstest_statistic" = ks_results$statistic[[1]]))



  }


  # Add country name and subject count into the results
  site_pvalues <- site_pvalues %>%
    left_join(y=site_metadata, by="site") %>%
    mutate(
      pvalue_kstest = ifelse(
        is.na(.data$pvalue_kstest) & as.numeric(.data$kstest_statistic) == 1,
        1e-100,
        .data$pvalue_kstest
      )
    )

  return(site_pvalues)

}


#' parse_readable_timeseries_combo_string
#'
#' Parse timepoint codes into a more human-readable format.
#'
#' @param timeseries_combo Semicolon-delimited string with time point codes.
#' @param this_parameter_id Parameter the time point combo is for.
#' @param this_timepoint_rank_to_name_mapping Mapping table from time point code to human-readable name.
#' @return Semicolon-delimited string of human-readable timepoint names.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
parse_readable_timeseries_combo_string <- function(timeseries_combo, this_parameter_id, this_timepoint_rank_to_name_mapping) {

  # Time point combination strings used by the script are not very human friendly.
  # This function converts these into more readable ones which one can show in the
  # TSO dashboard.

  output <- c()

  timepoints <- unlist( strsplit(timeseries_combo, ";") )

  for(this_timepoint in timepoints) {

    timepoint_components <- this_timepoint_rank_to_name_mapping %>%
      filter(.data$timepoint_rank == .env$this_timepoint & .data$parameter_id == .env$this_parameter_id) %>%
      pull(names)

    timepoint_components <- unlist( strsplit(timepoint_components, "_") )

    timepoint_components <- timepoint_components[timepoint_components != "ND"]

    timepoint_components <- paste(timepoint_components, collapse="_")

    output <- append(output, timepoint_components)

  }

  output <- as.character( paste(output, collapse=";") )

  return(output)

}

#' calculate_own_site_simil_score
#'
#' Calculates the own site similarity feature for each time series.
#' @inheritParams process_a_study
#' @param this_timeseries_wide Time series data in wide format.
#' @return Data frame with a feature value per subject.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
calculate_own_site_simil_score <- function(this_timeseries_wide, subjects) {

  this_timeseries_wide$subject_id <- as.character(this_timeseries_wide$subject_id)
  subj_keys <- this_timeseries_wide$subject_id
  this_timeseries_wide$subject_id <- NULL
  suppressWarnings (rownames(this_timeseries_wide) <- subj_keys)

  origvalues_distances <- dist(this_timeseries_wide)


  dist_object_labels <- attr(origvalues_distances, "Labels")
  dist_object_label_count <- attr(origvalues_distances, "Size")

  own_site_simil_scores <- rep(-1, times=dist_object_label_count)

  dist_object_label_sites <- data.frame("subject_id" = dist_object_labels) %>%
    left_join(y=select(subjects, c("subject_id", "site")), by="subject_id") %>%
    pull(.data$site)

  # Only process sites with more than one subject per time series
  sites_to_process <- subjects %>%
    filter(.data$site %in% dist_object_label_sites) %>%
    filter(.data$subject_id %in% dist_object_labels) %>%
    group_by(.data$site) %>%
    summarise(subj_count = n_distinct(.data$subject_id)) %>%
    filter(.data$subj_count > 1) %>%
    pull(.data$site)


  for(this_site in unique(sites_to_process)) {

    this_site_filter <- ifelse(dist_object_label_sites == this_site, 1, 0)

    site_subject_indices <- which(this_site_filter == 1)

    for(this_subject_index in site_subject_indices) {

      labels_before <- this_subject_index - 1
      labels_after <- dist_object_label_count - this_subject_index

      if(this_subject_index == 1) {

        dist_indices_before <- c()

        dist_indices_own_col <- seq(from=1, to=labels_after)

      } else if (this_subject_index == dist_object_label_count) {

        # Distances for the subject before the "own" subject's column.
        index_steps <- cumsum( seq(from=labels_after+labels_before-1, by=-1, length.out=labels_before-1) )
        dist_indices_before <- c(labels_before, labels_before+index_steps)

        dist_indices_own_col <- c()

      } else {

        # Distances for the subject before the "own" subject's column.
        index_steps <- cumsum( seq(from=labels_after+labels_before-1, by=-1, length.out=labels_before-1) )
        dist_indices_before <- c(labels_before, labels_before+index_steps)


        # Distances given in the subject's "own" column
        own_col_start_index <- sum(seq(from=dist_object_label_count-1, to=labels_after+1, by=-1)) + 1
        own_col_end_index <- own_col_start_index + labels_after - 1

        dist_indices_own_col <- seq(from=own_col_start_index, to=own_col_end_index)

      }


      distances_before <- origvalues_distances[dist_indices_before]
      distances_own_col <- origvalues_distances[dist_indices_own_col]

      distances <- c(distances_before, distances_own_col)

      distances_site_categ <- this_site_filter[-this_subject_index]

      own_site_simil_scores[this_subject_index] <- auroc(distances, distances_site_categ)

    }



  }


  same_site_simil_metric <- data.frame("subject_id" = dist_object_labels,
                                       "own_site_simil_score" = own_site_simil_scores) %>%
    filter(.data$own_site_simil_score >= 0)

  return(same_site_simil_metric)

}

#' auroc
#'
#' Fast function for calculating ROC AUC (original code by Miron Kursa https://mbq.me).
#'
#' @param score Score for ranking the observations.
#' @param bool Boolean giving the observation class.
#' @return ROC AUC value
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
auroc <- function(score, bool) {
  n1 <- sum(!bool)
  n2 <- sum(bool)
  U  <- sum(rank(-score)[!bool]) - n1 * (n1 + 1) / 2
  return(1 - U / n1 / n2)
}

#' calculate_ts_features
#'
#' Calculates features for a time series.
#'
#' @param this_timeseries_wide Time series data in wide format.
#' @param this_baseline Whether the data in this_timeseries_wide are original or baseline-adjusted measurements.
#' @param this_timeseries_features_to_calculate Vector of time series features to return.
#' @inheritParams process_a_study
#' @return Data frame with time series features.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
calculate_ts_features <- function(this_timeseries_wide, this_baseline, this_timeseries_features_to_calculate, subjects) {

  own_site_similarity_scores <- calculate_own_site_simil_score(this_timeseries_wide, subjects)

  subj_keys <- this_timeseries_wide$subject_id
  this_timeseries_wide$subject_id <- NULL

  timepoint_col_names <- colnames(this_timeseries_wide)

  averages <- apply(this_timeseries_wide, MARGIN = 1, mean, na.rm = TRUE)

  ts_features <- data.frame("subject_id" = subj_keys)

  if(ncol(this_timeseries_wide) > 1 & this_baseline == "original") {

    # Most time series features should only be calculated if more than time points are available

    standard_devs <- apply(this_timeseries_wide, MARGIN = 1, sd, na.rm = TRUE)
    unique_value_counts_relative <- apply(this_timeseries_wide, MARGIN = 1, function(x) n_distinct(x, na.rm = TRUE) / sum(!is.na(x)) )
    autocorrelations <- apply(this_timeseries_wide, MARGIN = 1, calculate_autocorrelation)
    value_ranges <- apply(this_timeseries_wide, MARGIN = 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE) )

    ts_features$sd <- standard_devs
    ts_features$average <- averages
    ts_features$autocorr <- autocorrelations
    ts_features$unique_value_count_relative <- unique_value_counts_relative
    ts_features$range <- value_ranges



  } else if (ncol(this_timeseries_wide) > 1 & this_baseline == "cfb") {

    # For change-from-baseline time series, return average and self-similarity to own site features.
    # Other features would be identical with the time series for original values.

    ts_features$average <- averages


  } else if (ncol(this_timeseries_wide) == 1) {

    # For a single time points, the only features that make sense are average (i.e. the only value itself) and
    # self-similarity to own site.

    ts_features$average <- averages

  }

  # Join the own site similarity scores and include only features that have been requested.
  ts_features <- ts_features %>%
    left_join(y=own_site_similarity_scores, by="subject_id") %>%
    pivot_longer(cols= - c("subject_id"), names_to="feature", values_to="value", values_drop_na = TRUE) %>%
    filter(.data$feature %in% this_timeseries_features_to_calculate)


  return(ts_features)

}


#' calculate_principal_components
#'
#' Calculate the first two principal components. These are used to visualize the relative distance between the time series/subjects.
#'
#' @param this_timeseries_wide Wide-format data frame with data to calculate the principal components for.
#' @return Data frame with the principal components.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
calculate_principal_components <- function(this_timeseries_wide) {

  subj_keys <- this_timeseries_wide$subject_id
  this_timeseries_wide$subject_id <- NULL

  # Impute missing values with row averages
  this_timeseries_wide <- as.matrix(this_timeseries_wide)
  k <- which(is.na(this_timeseries_wide), arr.ind=TRUE)
  this_timeseries_wide[k] <- rowMeans(this_timeseries_wide, na.rm=TRUE)[k[,1]]
  this_timeseries_wide <- as.data.frame(this_timeseries_wide)

  # First, remove constant columns since these will cause errors with prcomp
  this_timeseries_wide <- as.data.frame( this_timeseries_wide[,apply(this_timeseries_wide, 2, var, na.rm=TRUE) != 0] )

  principal_components <- data.frame()

  if(ncol(this_timeseries_wide) == 1) {
    # If we only have one time point, use the actual measurement as the first "principal component"
    # and assing the second to constant.

    principal_components <- this_timeseries_wide
    colnames(principal_components) <- "PC1"
    principal_components$PC2 <- 0

  }

  if(ncol(this_timeseries_wide) == 2) {
    # If we only have two time points, use the actual measurements as "principal components"

    principal_components <- this_timeseries_wide
    colnames(principal_components) <- c("PC1", "PC2")

  }

  if(ncol(this_timeseries_wide) > 2) {

    # For time series with more than two time points, represent the time series with the first
    # principal components.



    principal_components <- prcomp(this_timeseries_wide, scale. = FALSE)
    principal_components <- as.data.frame(  principal_components$x )
    principal_components <- principal_components[,1:2]

  }

  # Add subject IDs only if principal components were calculated
  if(ncol(principal_components) > 0) {

    principal_components$subject_id <- subj_keys

  }



  return(principal_components)

}


#' generate_wide_timeseries_table
#'
#' Transforms data for a time series as a wide data frame.
#'
#' @param this_parameter_id Parameter ID.
#' @param this_timepoints Semicolon-delimited string of time series timepoints.
#' @param this_subjects Semicolon-delimited string of subject IDs.
#' @param this_baseline If TRUE, the values should be baseline-adjusted.
#' @inheritParams process_a_study
#' @return Time series data in wide format.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
generate_wide_timeseries_table <- function(this_parameter_id, this_timepoints, this_subjects, this_baseline, data) {

  this_timepoints <- strsplit(this_timepoints, ";")[[1]]
  this_subjects <- strsplit(this_subjects, ";")[[1]]

  measurements <- data %>%
    ungroup() %>%
    filter(.data$timepoint_rank %in% this_timepoints & .data$parameter_id == this_parameter_id & .data$subject_id %in% this_subjects) %>%
    arrange(.data$timepoint_rank)

  # If the values should be change-from-baseline (cfb), remove the baseline from each result.
  if(this_baseline == "cfb") {

    measurements <- measurements %>%
      filter(!is.na(.data$baseline)) %>%
      mutate(result = .data$result - .data$baseline)

  }

  timeseries_wide <- measurements %>%
    pivot_wider(id_cols = c("subject_id"), names_from = "timepoint_rank", values_from = "result", values_fn = mean)

  return(timeseries_wide)

}

#' pick_subjects_for_custom_timeseries
#'
#' Picks subjects which have enough data for a time series.
#'
#' @param this_timepoints_and_subjects Data frame which tells which subjects have a result at which time point.
#' @param this_timepoints_string Time points that belong to the time series.
#' @param this_max_share_missing Maximum share of missing data points a subject can have.
#' @param this_parameter_id ID of the time series parameter.
#' @param this_baseline Whether the time series is for original or baseline-adjusted results.
#' @return Semicolon-delimited string of subject IDs.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
pick_subjects_for_custom_timeseries <- function(this_timepoints_and_subjects, this_timepoints_string, this_max_share_missing, this_parameter_id, this_baseline) {

  # If the time series should be baseline-adjusted, pick only subject-parameter pairs
  # with a baseline value available.
  if(this_baseline == "cfb") {

    this_timepoints_and_subjects <- this_timepoints_and_subjects %>%
      filter(.data$has_baseline_value == "Yes")

  }

  this_timepoints_string <- str_split(this_timepoints_string, pattern = ';')[[1]]

  timepoint_count <- length(this_timepoints_string)

  timepoint_combo_subjects <- this_timepoints_and_subjects %>%
    filter(.data$parameter_id == .env$this_parameter_id) %>%
    select(c("subject_id", "timepoint_rank")) %>%
    filter(.data$timepoint_rank %in% .env$this_timepoints_string) %>%
    group_by(.data$subject_id) %>%
    summarise(measurement_count = n()) %>%
    filter(.data$measurement_count >= ceiling( (1 - .env$this_max_share_missing) * .env$timepoint_count)) %>%
    pull(.data$subject_id)

  subject_ids_to_return <- paste(timepoint_combo_subjects, collapse=";")

}

#' pick_timepoint_combos
#'
#' Autogenerates one or more time series for the parameter.
#'
#' @param dataset Description of parameter x.
#' @param this_time_point_count_min Minimum number of time points a time series must have.
#' @param this_subject_count_min Minimum number of subjects per time series.
#' @param this_max_share_missing Maximum share of missing measuremnents a subject can have.
#' @param this_baseline Whether the time series is for actual measurements or change-from-baseline values.
#' @return Data frame with one or more time series for the parameter.
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
pick_timepoint_combos <- function(dataset, this_time_point_count_min, this_subject_count_min, this_max_share_missing, this_baseline) {


  timepoint_combos_to_return <- c()
  subject_ids_to_return <- c()

  prev_accepted_timeseries_subj_count <- 0

  # If the timepoint combos should be for baseline-adjusted values, only
  # include subject-parameter pairs for which a baseline value is available.
  if(this_baseline == "cfb") {

    dataset <- dataset %>%
      filter(.data$has_baseline_value == "Yes")

  }

  timepoint_ranks <- sort( unique(dataset$timepoint_rank) )

  # Execute this only if there are enough data points available.
  if(length(timepoint_ranks) >= this_time_point_count_min) {

    # Starting from last time point, define the time series as all time points from the first one to this_last_visit_index.
    # Keep time series which have enough eligible subjects and which are non-redundant.
    for(this_last_visit_index in length(timepoint_ranks):this_time_point_count_min) {

      this_ts_timepoints <- timepoint_ranks[1:this_last_visit_index]

      timepoint_count <- length(this_ts_timepoints)

      # List of subjects who have at most this_max_share_missing missing for the time points.
      timeseries_subjects <- dataset %>%
        filter(.data$timepoint_rank %in% .env$this_ts_timepoints) %>%
        group_by(.data$subject_id) %>%
        summarise(measurement_count = n()) %>%
        filter(.data$measurement_count >= ceiling( (1 - .env$this_max_share_missing) * .env$timepoint_count)) %>%
        pull(.data$subject_id)

      num_subjects_in_timeseries <- length(timeseries_subjects)

      if(num_subjects_in_timeseries >= this_subject_count_min) {

        # Only include time series with enough subjects and
        # if the number of subjects is at least 20 % greater than the subject count
        # for the previous time series.
        if(prev_accepted_timeseries_subj_count == 0 |
           ( num_subjects_in_timeseries - prev_accepted_timeseries_subj_count >= this_subject_count_min &
             num_subjects_in_timeseries / prev_accepted_timeseries_subj_count >= 1.2)
        ) {

          timepoint_combos_to_return <- append(timepoint_combos_to_return, paste(this_ts_timepoints, collapse=";"))
          subject_ids_to_return <- append(subject_ids_to_return, paste(timeseries_subjects, collapse=";"))

          prev_accepted_timeseries_subj_count <- num_subjects_in_timeseries

        }


      }



    }

  }

  return_df <- data.frame("timepoint_combo" = timepoint_combos_to_return, "timepoint_combo_subjects" = subject_ids_to_return)
  return(return_df)


}

#' check_input_data
#'
#' Makes sure the input data is in correct format etc. No return value is given. The function stops processing if something is not right.
#' @inheritParams process_a_study
#' @return no data returned
#'
#' @author Pekka Tiikkainen, \email{pekka.tiikkainen@@bayer.com}
check_input_data <- function(subjects, parameters, data, custom_timeseries, timeseries_features_to_calculate,
                             default_minimum_timepoints_per_series, default_minimum_subjects_per_series,
                             default_max_share_missing_timepoints_per_series, default_generate_change_from_baseline,
                             autogenerate_timeseries) {

  df_subjects_colnames_expected <- c("country", "subject_id", "site")

  df_parameters_colnames_expected <- c("parameter_id", "parameter_category_1", "parameter_category_2",
                                       "parameter_category_3", "parameter_name", "time_point_count_min", "subject_count_min", "max_share_missing",
                                       "generate_change_from_baseline")

  df_data_colnames_expected <- c("subject_id", "parameter_id", "timepoint_1_name", "timepoint_2_name", "timepoint_rank", "result", "baseline")

  df_custom_timeseries_colnamed_expected <- c("timeseries_id", "parameter_id", "timepoint_combo")

  allowed_timeseries_features <- c('autocorr', 'average', 'own_site_simil_score', 'sd', 'unique_value_count_relative', 'range')



  stopifnot("There is no data!" = nrow(data) > 0)
  stopifnot("There are no subjects!" = nrow(subjects) > 0)

  # Check that the input data frames have all the expected columns.
  stopifnot("The argument df 'subjects' does not have all the expected columns!" = setequal(colnames(subjects), df_subjects_colnames_expected))
  stopifnot("The argument df 'parameters' does not have all the expected columns!" = setequal(colnames(parameters), df_parameters_colnames_expected))
  stopifnot("The argument df 'data' does not have all the expected columns!" = setequal(colnames(data), df_data_colnames_expected))

  # And make sure those columns have correct data types. The is_logical test is for columns which
  # can be left empty and if there are no values, R labels these as "logical".

  # Data frame 'Subjects'
  stopifnot("Column 'subject_id' in 'subjects' must be a character!" = is.character(subjects$subject_id))
  stopifnot("Column 'country' in 'subjects' must be a character!" = is.character(subjects$country))
  stopifnot("Column 'site' in 'subjects' must be a character!" = is.character(subjects$site))

  # Data frame 'Parameters'
  stopifnot("Column 'parameter_id' in 'parameters' must be a character!" = is.character(parameters$parameter_id))
  stopifnot("Column 'parameter_category_1' in 'parameters' must be a character!" = is.character(parameters$parameter_category_1))
  stopifnot("Column 'parameter_category_2' in 'parameters' must be a character!" = is.character(parameters$parameter_category_2) | is.logical(parameters$parameter_category_2))
  stopifnot("Column 'parameter_category_3' in 'parameters' must be a character!" = is.character(parameters$parameter_category_3) | is.logical(parameters$parameter_category_3))
  stopifnot("Column 'parameter_name' in 'parameters' must be a character!" = is.character(parameters$parameter_name))
  stopifnot("Column 'time_point_count_min' in 'parameters' must be numeric!" = is.numeric(parameters$time_point_count_min) | is.logical(parameters$time_point_count_min))
  stopifnot("Column 'subject_count_min' in 'parameters' must be numeric!" = is.numeric(parameters$subject_count_min) | is.logical(parameters$subject_count_min))
  stopifnot("Column 'max_share_missing' in 'parameters' must be numeric!" = is.numeric(parameters$max_share_missing) | is.logical(parameters$max_share_missing))
  stopifnot("Column 'generate_change_from_baseline' in 'parameters' must be logical!" = is.logical(parameters$generate_change_from_baseline))

  # Data frame 'Data'
  stopifnot("Column 'subject_id' in 'data' must be a character!" = is.character(data$subject_id))
  stopifnot("Column 'parameter_id' in 'data' must be a character!" = is.character(data$parameter_id))
  stopifnot("Column 'timepoint_1_name' in 'data' must be a character!" = is.character(data$timepoint_1_name))
  stopifnot("Column 'timepoint_2_name' in 'data' must be a character!" = is.character(data$timepoint_2_name) | is.logical(data$timepoint_2_name))
  stopifnot("Column 'timepoint_rank' in 'data' must be numeric!" = is.numeric(data$timepoint_rank))
  stopifnot("Column 'result' in 'data' must be numeric!" = is.numeric(data$result))
  stopifnot("Column 'baseline' in 'data' must be numeric!" = is.numeric(data$baseline) | is.logical(data$baseline))

  if(nrow(custom_timeseries) > 0) {

    stopifnot("The argument df 'custom_timeseries' does not have all the expected columns!" = setequal(colnames(custom_timeseries), df_custom_timeseries_colnamed_expected))

    # Data frame 'custom_timeseries'
    stopifnot("Column 'timeseries_id' in 'custom_timeseries' must be a character!" = is.character(custom_timeseries$timeseries_id))
    stopifnot("Column 'parameter_id' in 'custom_timeseries' must be a character!" = is.character(custom_timeseries$parameter_id))
    stopifnot("Column 'timepoint_combo' in 'custom_timeseries' must be a character!" = is.character(custom_timeseries$timepoint_combo))

  }


  # Make sure that all feature names to calculate are correct.
  stopifnot("The argument 'timeseries_features_to_calculate' contains illegal values!" = all(timeseries_features_to_calculate %in% allowed_timeseries_features))

  # Make sure scalar arguments have correct data types.
  stopifnot("Argument default_minimum_timepoints_per_series is not numeric!" = is.numeric(default_minimum_timepoints_per_series))
  stopifnot("Argument default_minimum_subjects_per_series is not numeric!" = is.numeric(default_minimum_subjects_per_series))
  stopifnot("Argument default_max_share_missing_timepoints_per_series must be between 0 and 1!" = all(!is.na(default_max_share_missing_timepoints_per_series) & default_max_share_missing_timepoints_per_series >= 0 & default_max_share_missing_timepoints_per_series <= 1))
  stopifnot("Argument default_generate_change_from_baseline must be TRUE or FALSE!" = is.logical(default_generate_change_from_baseline))
  stopifnot("Argument autogenerate_timeseries must be TRUE or FALSE!" = is.logical(autogenerate_timeseries))

  # If the time series should not be auto generated, expect at least one custom timeseries.
  stopifnot("Custom timeseries must be defined if autogenerate_timeseries is set to FALSE!" = ifelse(autogenerate_timeseries, TRUE, nrow(custom_timeseries) > 0))

  # Stop if there are replicate subject ids in subjects
  replicate_subject_ids <- subjects %>%
    group_by(.data$subject_id) %>%
    summarise(rowcount = n()) %>%
    filter(.data$rowcount > 1) %>%
    pull(.data$subject_id)

  stopifnot("There are replicate subject IDs in the subjects df!" = length(replicate_subject_ids) == 0)

}


#' calculate_autocorrelation
#'
#' Calculates the autocorrelation feature for a time series. Default lag is one time point.
#'
#' @param this_timeseries Measurements of a time series for an individual subject ranked by time.
#' @param lag Lag parameter for auto-correlation. Default lag is one.
#' @return Autocorrelation coefficient.
calculate_autocorrelation <- function(this_timeseries, lag=1) {

  timeseries_length <- length(this_timeseries)

  auto_corr_coeff <- cor(this_timeseries[1:(timeseries_length-lag)],
                         this_timeseries[(1+lag):timeseries_length],
                         use = "pairwise.complete.obs")

  return(auto_corr_coeff)

}
