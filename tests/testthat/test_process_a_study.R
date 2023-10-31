test_that("process_a_study", {

  data("tsoa_data", package = "tsoa")

  ls_tsoa <- process_a_study(
    data = tsoa_data$data,
    subjects = tsoa_data$subjects,
    parameters = tsoa_data$parameters,
    custom_timeseries = tsoa_data$custom_timeseries,
    timeseries_features_to_calculate = c(
      "autocorr",
      "average",
      "own_site_simil_score",
      "sd",
      "unique_value_count_relative"
    ),
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

  expect_true(all(expected_tsoa_names %in%  names(ls_tsoa)))

  allow_na <- c(
    "fdr_corrected_pvalue_logp",
    "pvalue_kstest_logp"
  )

  expect_true(all(unlist(
    purrr::map(
      ls_tsoa,
      ~ ! any(is.na(
        select(., - any_of(allow_na))
        ))
      )
  )))

  expect_true(all(unlist(
    purrr::map(
      ls_tsoa,
      ~ nrow(.) > 0
    )
  )))

})


