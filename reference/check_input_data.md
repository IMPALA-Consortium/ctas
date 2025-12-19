# check_input_data

Makes sure the input data is in correct format etc. No return value is
given. The function stops processing if something is not right.

## Usage

``` r
check_input_data(
  subjects,
  parameters,
  data,
  custom_timeseries,
  custom_reference_groups,
  default_timeseries_features_to_calculate,
  default_minimum_timepoints_per_series,
  default_minimum_subjects_per_series,
  default_max_share_missing_timepoints_per_series,
  default_generate_change_from_baseline,
  autogenerate_timeseries
)
```

## Arguments

- subjects:

  Data frame with one row per study subject.

- parameters:

  Data frame with one row per study parameter (e.g. labs and vital
  signs) to process.

- data:

  Measurements of study parameters. One row per data point.

- custom_timeseries:

  Pre-defined time series.

- custom_reference_groups:

  Parameter-feature combinations for which the sites should be compared
  within the country or region - and not globally.

- default_timeseries_features_to_calculate:

  Semicolon-delimited string of time series features to calculate.

- default_minimum_timepoints_per_series:

  Default minimum timepoints per time series. This is used if no
  parameter-specific minimum is defined in parameters.

- default_minimum_subjects_per_series:

  Default minimum of subjects per time series. This is used if no
  parameter-specific minimum is defined in parameters.

- default_max_share_missing_timepoints_per_series:

  Maximum share of time points which can be missing per subject. This is
  used if no parameter-specific minimum is defined in parameters.

- default_generate_change_from_baseline:

  If set TRUE, time series and their features are calculate also for
  change-from-baseline values.

- autogenerate_timeseries:

  If set to TRUE, automatic definition of time series is used. If set to
  FALSE, custom_timeseries must have at least one time series defined.

## Value

no data returned

## Author

Pekka Tiikkainen, <pekka.tiikkainen@bayer.com>
