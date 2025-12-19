# pick_subjects_for_custom_timeseries

Picks subjects which have enough data for a time series.

## Usage

``` r
pick_subjects_for_custom_timeseries(
  this_timepoints_and_subjects,
  this_timepoints_string,
  this_max_share_missing,
  this_parameter_id,
  this_baseline
)
```

## Arguments

- this_timepoints_and_subjects:

  Data frame which tells which subjects have a result at which time
  point.

- this_timepoints_string:

  Time points that belong to the time series.

- this_max_share_missing:

  Maximum share of missing data points a subject can have.

- this_parameter_id:

  ID of the time series parameter.

- this_baseline:

  Whether the time series is for original or baseline-adjusted results.

## Value

Semicolon-delimited string of subject IDs.

## Author

Pekka Tiikkainen, <pekka.tiikkainen@bayer.com>
