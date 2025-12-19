# calculate_ts_features

Calculates features for a time series.

## Usage

``` r
calculate_ts_features(
  this_timeseries_wide,
  this_baseline,
  this_timeseries_features_to_calculate,
  subjects
)
```

## Arguments

- this_timeseries_wide:

  Time series data in wide format.

- this_baseline:

  Whether the data in this_timeseries_wide are original or
  baseline-adjusted measurements.

- this_timeseries_features_to_calculate:

  Semicolon-delimited list of time series features to return.

- subjects:

  Data frame with one row per study subject.

## Value

Data frame with time series features.

## Author

Pekka Tiikkainen, <pekka.tiikkainen@bayer.com>
