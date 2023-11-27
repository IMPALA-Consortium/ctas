test_that("calculate_autocorrelation tolerates NA", {
  set.seed(1)
  expect_true(! is.na(calculate_autocorrelation(c(NA,
                                                  rnorm(10, 5, 1),
                                                  NA, rnorm(10, 5, 1),
                                                  NA))))
})


test_that("ks.test not to return NA", {
  this_data <- tibble(
    site = c("A", "B"),
    country = c("C", "C"),
    region = c("D", "D"),
    value = c(list(rnorm(1000, 5, 0.1)), list(rnorm(1000, 50, 0.1)))
  ) %>%
    unnest(value) %>%
    mutate(
      subject_id = as.character(row_number())
    )
  
  this_feature = "test"
  this_ref_group = "global"
  
  df_kstest <- calculate_site_bias_ts_features(this_feature, this_data, this_ref_group)

})
