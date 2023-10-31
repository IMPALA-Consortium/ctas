test_that("calculate_autocorrelation tolerates NA", {
  set.seed(1)
  expect_true(! is.na(calculate_autocorrelation(c(NA,
                                                  rnorm(10, 5, 1),
                                                  NA, rnorm(10, 5, 1),
                                                  NA))))
})
