
feats <- c(
  "autocorr",
  "average",
  "own_site_simil_score",
  "sd",
  "unique_value_count_relative",
  "lof",
  "range"
)

feats_collapse = paste(feats, collapse = ";")

test_that("process_a_study", {

  data("tsoa_data", package = "tsoa")

  ls_tsoa <- process_a_study(
    data = tsoa_data$data,
    subjects = tsoa_data$subjects,
    parameters = tsoa_data$parameters,
    custom_timeseries = tsoa_data$custom_timeseries,
    custom_reference_groups = tsoa_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 3,
    default_minimum_subjects_per_series = 3,
    default_max_share_missing_timepoints_per_series = 0.5,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE
  )

  expected_tsoa_names <- c(
    "site_scores",
    "PCA_coordinates",
    "timeseries_features",
    "timeseries"
  )

  expect_true(all(feats %in% unique(ls_tsoa$site_scores$feature)))

  expect_true(all(expected_tsoa_names %in%  names(ls_tsoa)))

  expect_true(all(unlist(
    purrr::map(
      ls_tsoa,
      ~ ! any(is.na(.)
        ))
      )
  ))

  expect_true(all(unlist(
    purrr::map(
      ls_tsoa,
      ~ nrow(.) > 0
    )
  )))

})

test_that("process_a_study - optimize_sites_and_patients, backwards compatibility", {

  data("tsoa_data", package = "tsoa")

  # only keep the first 25 percent of timepoint ranks for site AAA
  tsoa_data$data <- tsoa_data$data %>%
    left_join(
      tsoa_data$subjects %>%
        select(subject_id, site),
      by = "subject_id"
    ) %>%
    filter(
      site != "AAA" | percent_rank(timepoint_rank) <= 0.25,
      .by = "site"
    ) %>%
    select(-site)

  # in default we expect site AAA which has only patients with
  # a low number of patients to get dropped from the analysis
  ls_tsoa_default <- process_a_study(
    data = tsoa_data$data,
    subjects = tsoa_data$subjects,
    parameters = tsoa_data$parameters,
    custom_timeseries = tsoa_data$custom_timeseries,
    custom_reference_groups = tsoa_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 3,
    default_minimum_subjects_per_series = 3,
    default_max_share_missing_timepoints_per_series = 0.5,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE
  )

  # if optimize_sites_and_patients an additional timeseries will be created
  # that includes the maximum number of sites and patients thus keeping site AAA
  ls_tsoa_opt <- process_a_study(
    data = tsoa_data$data,
    subjects = tsoa_data$subjects,
    parameters = tsoa_data$parameters,
    custom_timeseries = tsoa_data$custom_timeseries,
    custom_reference_groups = tsoa_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 3,
    default_minimum_subjects_per_series = 3,
    default_max_share_missing_timepoints_per_series = 0.5,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE,
    optimize_sites_and_patients = TRUE
  )

  expect_true(! "AAA" %in% ls_tsoa_default$site_scores$site)
  expect_true("AAA" %in% ls_tsoa_opt$site_scores$site)
  expect_true(all(ls_tsoa_default$timeseries$timepoint_combo %in% ls_tsoa_opt$timeseries$timepoint_combo))
  expect_true(! all(ls_tsoa_opt$timeseries$timepoint_combo %in% ls_tsoa_default$timeseries$timepoint_combo))

})




