# process_a_study

Main function for calculating time series features for study parameters,
and flagging study sites with systematic bias.

## Usage

``` r
process_a_study(
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
  autogenerate_timeseries,
  optimize_sites_and_patients = FALSE,
  site_scoring_method = "ks"
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

- optimize_sites_and_patients:

  If set to TRUE, always creates timeseries with as many sites and
  patients as possible while respecting the other function parameters.
  Default:FALSE

- site_scoring_method:

  How to score sites ("ks" = Kolmogorov-Smirnov, "mixedeffects" = mixed
  effects modelling, "avg_feat_value" = Average site feature value.
  Default:ks

## Value

List with four data frames. Timeseries: definition of the time series
used. Timeseries_features: features calculated from the time series.
PCA_coordinates: principal components of individual time series for
visualizing similarity. Site_scores: biasness scores for sites.

## Author

Pekka Tiikkainen, <pekka.tiikkainen@bayer.com>
