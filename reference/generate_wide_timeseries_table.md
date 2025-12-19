# generate_wide_timeseries_table

Transforms data for a time series as a wide data frame.

## Usage

``` r
generate_wide_timeseries_table(
  this_parameter_id,
  this_timepoints,
  this_subjects,
  this_baseline,
  data
)
```

## Arguments

- this_parameter_id:

  Parameter ID.

- this_timepoints:

  Semicolon-delimited string of time series timepoints.

- this_subjects:

  Semicolon-delimited string of subject IDs.

- this_baseline:

  If TRUE, the values should be baseline-adjusted.

- data:

  Measurements of study parameters. One row per data point.

## Value

Time series data in wide format.

## Author

Pekka Tiikkainen, <pekka.tiikkainen@bayer.com>
