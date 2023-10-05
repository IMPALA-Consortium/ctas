## code to prepare `DATASET` dataset goes here

library(dplyr)
library(tibble)
library(purrr)
library(tidyr)

timepoint_names <- combn(LETTERS, 2, FUN = paste, collapse = "")


ts <- tibble::tibble(
    site = as.character(seq(1, 35))
  ) %>%
  dplyr::mutate(
    country = sample(seq(1, 5), size = nrow(.), replace = TRUE),
    country = LETTERS[country],
    subject_id = rpois(nrow(.), lambda = 5),
    subject_id = purrr::map(subject_id, ~ seq(1, .))
  ) %>%
  tidyr::unnest(subject_id)  %>%
  dplyr::mutate(
    subject_id = row_number(),
    subject_id = as.character(subject_id)
  ) %>%
  dplyr::mutate(
    timepoint_rank = rpois(nrow(.), lambda = 20),
    timepoint_rank = purrr::map(timepoint_rank, ~ seq(1, .))
  ) %>%
  tidyr::unnest(timepoint_rank) %>%
  dplyr::mutate(
    timepoint_1_name = timepoint_names[timepoint_rank],
  )

get_parameter <- function(ts, name) {

  rnorm_param <- ts %>%
    distinct(subject_id) %>%
    mutate(
      avg = rnorm(nrow(.), 30, 5),
      sd = runif(nrow(.), min = 1, max = 10),
    )

  ts %>%
    left_join(rnorm_param, by = "subject_id") %>%
    mutate(
      result = purrr::map2_dbl(avg, sd, function(x,y) rnorm(1, x, y)),
      result = ifelse(runif(nrow(.), 0, 1) <= 0.3, NA, result),
      parameter_id = name
    ) %>%
    select(- avg, - sd)
}

param1 <- get_parameter(ts, "param1")
param2 <- get_parameter(ts, "param2")

data <- bind_rows(
    param1,
    param2
  ) %>%
  mutate(
    timepoint_2_name = "timepoint name 2",
    baseline = NA
  ) %>%
  select(- site, - country)


parameters <- data %>%
  distinct(parameter_id) %>%
  mutate(
    parameter_name = parameter_id,
    parameter_category_1 = "category 1",
    parameter_category_2 = "category 2",
    parameter_category_3 = "category 3",
    time_point_count_min = NA,
    subject_count_min = NA,
    max_share_missing = NA,
    generate_change_from_baseline = NA
  )

subjects <- ts %>%
  distinct(subject_id, site, country)

custom_timeseries <- tibble(
  timeseries_id = NULL,
  parameter_id = NULL,
  timepoint_combo = NULL
)

tsoa_data <- list(
  data = data,
  parameters = parameters,
  subjects = subjects,
  custom_timeseries = custom_timeseries
)

usethis::use_data(tsoa_data, overwrite = TRUE)
