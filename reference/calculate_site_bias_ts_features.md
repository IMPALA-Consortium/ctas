# calculate_site_bias_ts_features

Performs Kolmogorov-Smirnov (KS) test for a particular parameter and
time series feature.

## Usage

``` r
calculate_site_bias_ts_features(this_feature, this_data, this_ref_group)
```

## Arguments

- this_feature:

  Type of the time series feature.

- this_data:

  Feature values for each subject.

- this_ref_group:

  Reference group for the sites (country, region or global).

## Value

Data frame with KS p-values for each site.

## Author

Pekka Tiikkainen, <pekka.tiikkainen@bayer.com>
