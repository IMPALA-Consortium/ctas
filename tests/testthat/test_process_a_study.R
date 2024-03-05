
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

  data("ctas_data", package = "ctas")

  ls_ctas <- process_a_study(
    data = ctas_data$data,
    subjects = ctas_data$subjects,
    parameters = ctas_data$parameters,
    custom_timeseries = ctas_data$custom_timeseries,
    custom_reference_groups = ctas_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 3,
    default_minimum_subjects_per_series = 3,
    default_max_share_missing_timepoints_per_series = 0.5,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE
  )

  expected_ctas_names <- c(
    "site_scores",
    "PCA_coordinates",
    "timeseries_features",
    "timeseries"
  )

  expect_true(all(feats %in% unique(ls_ctas$site_scores$feature)))

  expect_true(all(expected_ctas_names %in%  names(ls_ctas)))

  expect_true(all(unlist(
    purrr::map(
      ls_ctas,
      ~ ! any(is.na(.)
        ))
      )
  ))

  expect_true(all(unlist(
    purrr::map(
      ls_ctas,
      ~ nrow(.) > 0
    )
  )))

})

test_that("process_a_study - optimize_sites_and_patients, backwards compatibility", {

  data("ctas_data", package = "ctas")

  # only keep the first 25 percent of timepoint ranks for site AAA
  ctas_data$data <- ctas_data$data %>%
    left_join(
      ctas_data$subjects %>%
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
  ls_ctas_default <- process_a_study(
    data = ctas_data$data,
    subjects = ctas_data$subjects,
    parameters = ctas_data$parameters,
    custom_timeseries = ctas_data$custom_timeseries,
    custom_reference_groups = ctas_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 3,
    default_minimum_subjects_per_series = 3,
    default_max_share_missing_timepoints_per_series = 0.5,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE
  )

  # if optimize_sites_and_patients an additional timeseries will be created
  # that includes the maximum number of sites and patients thus keeping site AAA
  expect_warning(
      ls_ctas_opt <- process_a_study(
      data = ctas_data$data,
      subjects = ctas_data$subjects,
      parameters = ctas_data$parameters,
      custom_timeseries = ctas_data$custom_timeseries,
      custom_reference_groups = ctas_data$custom_reference_groups,
      default_timeseries_features_to_calculate = feats_collapse,
      default_minimum_timepoints_per_series = 3,
      default_minimum_subjects_per_series = 3,
      default_max_share_missing_timepoints_per_series = 0.5,
      default_generate_change_from_baseline = FALSE,
      autogenerate_timeseries = TRUE,
      optimize_sites_and_patients = TRUE
    ),
    regexp = "NA in dist object returning NA for lof"
  )

  expect_true(! "AAA" %in% ls_ctas_default$site_scores$site)
  expect_true("AAA" %in% ls_ctas_opt$site_scores$site)
  expect_true(all(ls_ctas_default$timeseries$timepoint_combo %in% ls_ctas_opt$timeseries$timepoint_combo))
  expect_true(! all(ls_ctas_opt$timeseries$timepoint_combo %in% ls_ctas_default$timeseries$timepoint_combo))

})


test_that("process_a_study - default_max_share_missing_timepoints_per_series <- 0", {

  data("ctas_data", package = "ctas")

  ctas_data$data <- ctas_data$data %>%
    mutate(
      result = ifelse(is.na(result), mean(result, na.rm = TRUE), result),
      # only timepoints that are not null for at least one patient are considered
      result = ifelse(timepoint_rank == 1 & subject_id != "1", NA, result)
    )

  ls_ctas <- process_a_study(
    data = ctas_data$data,
    subjects = ctas_data$subjects,
    parameters = ctas_data$parameters,
    custom_timeseries = ctas_data$custom_timeseries,
    custom_reference_groups = ctas_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 1,
    default_minimum_subjects_per_series = 3,
    default_max_share_missing_timepoints_per_series = 0,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE
  )

  expect_true(all(purrr::map_lgl(ls_ctas, is.null)))

})

test_that("process_a_study - default_minimum_timepoints_per_series <- 1e6", {

  data("ctas_data", package = "ctas")

  ls_ctas <- process_a_study(
    data = ctas_data$data,
    subjects = ctas_data$subjects,
    parameters = ctas_data$parameters,
    custom_timeseries = ctas_data$custom_timeseries,
    custom_reference_groups = ctas_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 1e6,
    default_minimum_subjects_per_series = 3,
    default_max_share_missing_timepoints_per_series = 1,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE
  )

  expect_true(all(purrr::map_lgl(ls_ctas, is.null)))

})

test_that("process_a_study - default_minimum_subjects_per_series <- 1e6", {

  data("ctas_data", package = "ctas")

  ls_ctas <- process_a_study(
    data = ctas_data$data,
    subjects = ctas_data$subjects,
    parameters = ctas_data$parameters,
    custom_timeseries = ctas_data$custom_timeseries,
    custom_reference_groups = ctas_data$custom_reference_groups,
    default_timeseries_features_to_calculate = feats_collapse,
    default_minimum_timepoints_per_series = 1,
    default_minimum_subjects_per_series = 1e6,
    default_max_share_missing_timepoints_per_series = 1,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE
  )

  expect_true(all(purrr::map_lgl(ls_ctas, is.null)))

})

test_that("process_a_study - default_minimum_subjects_per_series <- 1", {

  data("ctas_data", package = "ctas")

  df_lonely <- ctas_data$data %>%
    filter(subject_id == "1", parameter_id == "param1") %>%
    mutate(timepoint_rank = timepoint_rank + max(timepoint_rank))

  ctas_data$data <- bind_rows(
    ctas_data$data,
    df_lonely
  )

  expect_error(
    ls_ctas <- process_a_study(
      data = ctas_data$data,
      subjects = ctas_data$subjects,
      parameters = ctas_data$parameters,
      custom_timeseries = ctas_data$custom_timeseries,
      custom_reference_groups = ctas_data$custom_reference_groups,
      default_timeseries_features_to_calculate = feats_collapse,
      default_minimum_timepoints_per_series = 1,
      default_minimum_subjects_per_series = 1,
      default_max_share_missing_timepoints_per_series = 1,
      default_generate_change_from_baseline = FALSE,
      autogenerate_timeseries = TRUE
    ),
    "Minimum value for default_minimum_subjects_per_series is two!"
  )

})

test_that("process_a_study - avoid lof() error - minPts has to be at least 2 and not larger than the number of points", {

  data <- tibble(
    STUDY = rep("A", 372),
    PATNUM = c(rep("E", 10), rep("F", 23), rep("A", 36), rep("B", 80), rep("C", 158), rep("D", 29), rep("E", 36)),
    SITE_NUMBER = c(rep("D", 10), rep("E", 23), rep("C", 36), rep("A", 80), rep("B", 158), rep("E", 29), rep("D", 36)),
    TIMEPOINT_NAME = rep("A", 372),
    TIMEPOINT_RANK = c(11:21, 1:23, 1:36, 1:80, 1:158, 1:29, 1:35),
    NAME = rep("A", 372),
    PARAMETER_CATEGORY_1 = rep("A", 372),
    PARAMETER_CATEGORY_2 = rep("A", 372),
    PARAMETER_CATEGORY_3 = rep("A", 372),
    VALUE = sample(c(1:1000, NA), 372, replace = TRUE),
    COUNTRY = c(rep("B", 336), rep("A", 10), rep("E", 26)),
    REGION = rep("B", 372),
    SPLIT_BY = rep("A", 372)
  ) %>%
    rename_with(tolower) %>%
    rename(
      parameter_id = split_by,
      parameter_name = name,
      site = site_number,
      subject_id = patnum,
      timepoint_1_name = timepoint_name,
      result = value
    ) %>%
    filter(! (subject_id == "D" & timepoint_rank > 50))

  df_parameters <- data %>%
    distinct(
      .data$parameter_id,
      .data$parameter_name,
      .data$parameter_category_1,
      .data$parameter_category_2,
      .data$parameter_category_3
    ) %>%
    mutate(
      timeseries_features_to_calculate = NA,
      use_only_custom_timeseries = FALSE,
      time_point_count_min = NA,
      subject_count_min = NA,
      max_share_missing = NA,
      generate_change_from_baseline = NA
    )
  # the error occurs for sites with only one patient
  df_subjects <- data %>%
    count(.data$subject_id, .data$country, .data$region, .data$site) %>%
    arrange(.data$subject_id, desc(.data$n)) %>%
    filter(row_number() == 1, .by = "subject_id") %>%
    select(- n)

  df_data <- data %>%
    distinct(
      .data$subject_id,
      .data$parameter_id,
      .data$timepoint_1_name,
      .data$timepoint_rank,
      .data$result
    ) %>%
    mutate(
      timepoint_2_name = NA,
      baseline = NA
    )

  df_custom_timeseries <- tibble(
    timeseries_id = character(),
    parameter_id = character(),
    timepoint_combo = character()
  )

  df_custom_reference_groups <- tibble(
    parameter_id = character(),
    feature = character(),
    ref_group = character()
  )

  ls_ctas <- ctas::process_a_study(
    subjects = df_subjects,
    parameters = df_parameters,
    data = df_data,
    custom_timeseries = df_custom_timeseries,
    custom_reference_groups = df_custom_reference_groups,
    default_timeseries_features_to_calculate = "lof",
    default_minimum_timepoints_per_series = 1,
    default_minimum_subjects_per_series = 2,
    default_max_share_missing_timepoints_per_series = 0.4,
    default_generate_change_from_baseline = FALSE,
    autogenerate_timeseries = TRUE,
    optimize_sites_and_patients = TRUE
  )

  # data has very few subjects we expect score to be 0
  expect_true(all(ls_ctas$site_scores$fdr_corrected_pvalue_logp == 0))

})
