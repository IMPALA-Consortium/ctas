<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->
<!-- badges: end -->

# Time Series Outliers and Anomalies (TSOA)

Pekka Tiikkainen, Bayer Pharmaceuticals

<!-- badges: start -->
<!-- badges: end -->

## Introduction

TSOA is an R package for identifying outliers and anomalies in clinical
trial time series. Its main focus is on flagging sites with one or more
study parameters whose profiles differ from those of the other sites.
However, also anomalies for individual subjects can be identified with
the results the package gives.

This document focuses on the data format requirements of the input and
output data. For details on methodology, please refer to the paper and
presentation published at the PHUSE EU Connect 22 meeting. The files are
attached to this repository.

## How to use the package

First the clinical data has to be wrangled in a format compatible with
the main function’s parameters (please see below for definitions).
Please note that data preparation step is intentionally out-of-scope for
the package because each user is likely to have their data in a
different format.

Then the main function (process_a_study) is called separately for each
study you wish to process. The function returns the results as a list of
data frames. For details on the output, please see below.

## Example

``` r

library(tsoa)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

data("tsoa_data", package = "tsoa")

tsoa_data
#> $data
#> # A tibble: 13,774 × 7
#>    subject_id timepoint_rank timepoint_1_name result parameter_id timepoint_2_name baseline
#>    <chr>               <int> <chr[1d]>         <dbl> <chr>        <chr>            <lgl>   
#>  1 1                       1 AB                 35.5 param1       timepoint name 2 NA      
#>  2 1                       2 AC                 NA   param1       timepoint name 2 NA      
#>  3 1                       3 AD                 29.2 param1       timepoint name 2 NA      
#>  4 1                       4 AE                 NA   param1       timepoint name 2 NA      
#>  5 1                       5 AF                 NA   param1       timepoint name 2 NA      
#>  6 1                       6 AG                 31.2 param1       timepoint name 2 NA      
#>  7 1                       7 AH                 33.7 param1       timepoint name 2 NA      
#>  8 1                       8 AI                 27.5 param1       timepoint name 2 NA      
#>  9 1                       9 AJ                 NA   param1       timepoint name 2 NA      
#> 10 1                      10 AK                 32.9 param1       timepoint name 2 NA      
#> # … with 13,764 more rows
#> 
#> $parameters
#> # A tibble: 2 × 11
#>   parameter_id parameter_name parameter_category_1 parameter_category_2 parameter_category_3 time_po…¹ subje…² max_s…³ gener…⁴ times…⁵ use_o…⁶
#>   <chr>        <chr>          <chr>                <chr>                <chr>                <lgl>     <lgl>   <lgl>   <lgl>   <lgl>   <lgl>  
#> 1 param1       param1         category 1           category 2           category 3           NA        NA      NA      NA      NA      FALSE  
#> 2 param2       param2         category 1           category 2           category 3           NA        NA      NA      NA      NA      FALSE  
#> # … with abbreviated variable names ¹​time_point_count_min, ²​subject_count_min, ³​max_share_missing, ⁴​generate_change_from_baseline,
#> #   ⁵​timeseries_features_to_calculate, ⁶​use_only_custom_timeseries
#> 
#> $subjects
#> # A tibble: 349 × 4
#>    subject_id site  country region
#>    <chr>      <chr> <chr>   <chr> 
#>  1 1          AAA   AA      A     
#>  2 2          AAB   AA      A     
#>  3 3          AAB   AA      A     
#>  4 4          AAB   AA      A     
#>  5 5          AAB   AA      A     
#>  6 6          AAC   AA      A     
#>  7 7          AAC   AA      A     
#>  8 8          AAC   AA      A     
#>  9 9          AAC   AA      A     
#> 10 10         AAC   AA      A     
#> # … with 339 more rows
#> 
#> $custom_timeseries
#> # A tibble: 0 × 0
#> 
#> $custom_reference_groups
#> # A tibble: 0 × 0


ls_tsoa <- process_a_study(
  data = tsoa_data$data,
  subjects = tsoa_data$subjects,
  parameters = tsoa_data$parameters,
  custom_timeseries = tsoa_data$custom_timeseries,
  custom_reference_groups = tsoa_data$custom_reference_groups,
  default_timeseries_features_to_calculate = c(
    "autocorr",
    "average",
    "own_site_simil_score",
    "sd",
    "unique_value_count_relative",
    "lof",
    "range"
  ),
  default_minimum_timepoints_per_series = 3,
  default_minimum_subjects_per_series = 3,
  default_max_share_missing_timepoints_per_series = 0.5,
  default_generate_change_from_baseline = FALSE,
  autogenerate_timeseries = TRUE
)

ls_tsoa
#> $timeseries
#> # A tibble: 10 × 6
#>    timeseries_id          parameter_id baseline timepoint_combo                                                                timep…¹ timep…²
#>    <chr>                  <chr>        <chr>    <chr>                                                                          <chr>     <dbl>
#>  1 ts_1_autogen_original  param1       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29… AB_tim…      32
#>  2 ts_2_autogen_original  param1       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28     AB_tim…      28
#>  3 ts_3_autogen_original  param1       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24                 AB_tim…      24
#>  4 ts_4_autogen_original  param1       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20                             AB_tim…      20
#>  5 ts_5_autogen_original  param1       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14                                               AB_tim…      14
#>  6 ts_6_autogen_original  param2       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29… AB_tim…      33
#>  7 ts_7_autogen_original  param2       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28;29… AB_tim…      30
#>  8 ts_8_autogen_original  param2       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24;25;26;27;28     AB_tim…      28
#>  9 ts_9_autogen_original  param2       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20;21;22;23;24                 AB_tim…      24
#> 10 ts_10_autogen_original param2       original 1;2;3;4;5;6;7;8;9;10;11;12;13;14;15;16;17;18;19;20                             AB_tim…      20
#> # … with abbreviated variable names ¹​timepoint_combo_readable, ²​timepoint_count
#> 
#> $timeseries_features
#> # A tibble: 292 × 7
#>    timeseries_id         subject_id feature  feature_value site  country region
#>    <chr>                 <chr>      <chr>            <dbl> <chr> <chr>   <chr> 
#>  1 ts_1_autogen_original 1          autocorr        0.354  AAA   AA      A     
#>  2 ts_1_autogen_original 2          autocorr       -0.397  AAA   AA      A     
#>  3 ts_1_autogen_original 3          autocorr        0.119  AAA   AA      A     
#>  4 ts_1_autogen_original 4          autocorr        0.389  AAA   AA      A     
#>  5 ts_1_autogen_original 17         autocorr        0.309  AAC   AA      A     
#>  6 ts_1_autogen_original 19         autocorr        0.142  BAA   BA      B     
#>  7 ts_1_autogen_original 20         autocorr        0.0376 BAA   BA      B     
#>  8 ts_1_autogen_original 21         autocorr       -0.289  BAA   BA      B     
#>  9 ts_1_autogen_original 25         autocorr       -0.317  BAB   BA      B     
#> 10 ts_1_autogen_original 31         autocorr        0.0678 BAC   BA      B     
#> # … with 282 more rows
#> # ℹ Use `print(n = ...)` to see more rows
#> 
#> $PCA_coordinates
#> # A tibble: 292 × 4
#>    timeseries_id         subject_id    pc1    pc2
#>    <chr>                 <chr>       <dbl>  <dbl>
#>  1 ts_1_autogen_original 1          -20.8    6.00
#>  2 ts_1_autogen_original 2          -33.1    7.02
#>  3 ts_1_autogen_original 3           21.4   -2.27
#>  4 ts_1_autogen_original 4           31.0   -2.94
#>  5 ts_1_autogen_original 17          32.8   -3.77
#>  6 ts_1_autogen_original 19          -9.13   2.13
#>  7 ts_1_autogen_original 20          16.0   -9.33
#>  8 ts_1_autogen_original 21         -15.2    3.50
#>  9 ts_1_autogen_original 25         -22.4  -14.6 
#> 10 ts_1_autogen_original 31          29.8   28.6 
#> # … with 282 more rows
#> # ℹ Use `print(n = ...)` to see more rows
#> 
#> $site_scores
#> # A tibble: 83 × 10
#>    timeseries_id         site  country region feature  pvalue_kstest_logp kstest_statistic fdr_corrected_pvalue_logp ref_group subject_count
#>    <chr>                 <chr> <chr>   <chr>  <chr>                 <dbl>            <dbl>                     <dbl> <chr>             <int>
#>  1 ts_1_autogen_original AAA   AA      A      autocorr            0.393              0.5                           0 global                4
#>  2 ts_1_autogen_original AAC   AA      A      autocorr            0.426              0.867                         0 global                1
#>  3 ts_1_autogen_original BAA   BA      B      autocorr            0.135              0.385                         0 global                3
#>  4 ts_1_autogen_original BAB   BA      B      autocorr            0.602              0.933                         0 global                1
#>  5 ts_1_autogen_original BAC   BA      B      autocorr            0.0580             0.6                           0 global                1
#>  6 ts_1_autogen_original CAA   CA      C      autocorr            0.00383            0.25                          0 global                4
#>  7 ts_1_autogen_original CAC   CA      C      autocorr            0.155              0.5                           0 global                2
#>  8 ts_2_autogen_original AAA   AA      A      autocorr            0.338              0.389                         0 global                6
#>  9 ts_2_autogen_original AAB   AA      A      autocorr            0.115              0.455                         0 global                2
#> 10 ts_2_autogen_original AAC   AA      A      autocorr            0.301              0.783                         0 global                1
#> # … with 73 more rows
#> # ℹ Use `print(n = ...)` to see more rows

```

## Input parameters

The function process_a_study has a number of parameters that need to be
provided. Some of these are data frames while others are simple
variables.

The following gives more details on each parameter. Examples are also
provided where needed.

### default_timeseries_features_to_calculate

Vector of features to calculate for each time series by default. 
The vector must contain at least one of the following feature codes:

| Feature code                | Description                                                      |
|------------------------------------|------------------------------------|
| autocorr                    | Auto-correlation                                                 |
| average                     | Average value                                                    |
| own_site_simil_score        | Measure of co-clustering of time series from the same study site |
| sd                          | Standard deviation                                               |
| unique_value_count_relative | Number of unique values divided by number of values available    |
| range                       | The range of values in a time series                             |
| lof                         | Local outlier factor                                             |

### default_minimum_timepoints_per_series

Minimum number of time points for auto-generated time series. This value
will be used for a parameter if they do not have the
time_point_count_min column defined in the parameters df.

### default_minimum_subjects_per_series

Minimum number of eligible subjects for auto-generated time series. This
value will be used for a parameter if they do not have the
subject_count_min column defined in the parameters df.

### default_max_share_missing_timepoints_per_series

Maximum share of missing time points a subject can have in order to be
eligible for a time series. This value will be used for a parameter if
they do not have the max_share_missing defined in the parameters df.

### default_generate_change_from_baseline

If TRUE, change-from-baseline adjusted version of each time series is
generated (in addition to the time series with actual measurements).

### autogenerate_timeseries

If TRUE, time series are autogenerated in a data-driven matter for all
parameters. If set to FALSE, the custom_timeseries data frame must have
at least one entry.

### subjects

Data frame with one row per study subject.

Data frame columns:

| Column name | Data type | Data required | Description                            |
|-------------|-----------|---------------|----------------------------------------|
| subject_id  | chr       | Y             | Unique identifier for each subject.     |
| country     | chr       | Y             | Country where the subject is enrolled. |
| region      | chr       | Y             | Region of the world the site belongs to.|
| site        | chr       | Y             | Name of the study site.                |

### parameters

Data frame with one row for each parameter that should be processed.
Parameters are, for example, laboratory assays and vital sign
measurements done during the course of a study.

Data frame columns:

| Column name                   | Data type | Data required | Description                                                                                                                                       |
|------------------|------------------|------------------|------------------|
| parameter_id                  | chr       | Y             | Parameter identifier                                                                                                                              |
| parameter_category_1          | chr       |               | Top level category of the parameter (e.g. LB for a laboratory assay).                                                                             |
| parameter_category_2          | chr       |               | Second level category of the parameter (e.g. General chemistry for a laboratory assay).                                                           |
| parameter_category_3          | chr       |               | Third level category of the parameter (e.g. Local lab for a laboratory assay).                                                                    |
| parameter_name                | chr       | Y             | Parameter name                                                                                                                                    |
| time_point_count_min          | int       |               | Minimum number of time points required for auto-generated time series. If not provided, global default value is used.                             |
| subject_count_min             | int       |               | Minimum number of eligible subjects required for auto-generated time series. If not provided, global default value is used.                       |
| max_share_missing             | float     |               | Maximum number of missing data points an eligible subject can have for auto-generated time series. If not provided, global default value is used. |
| generate_change_from_baseline | boolean   |               | Whether to also generate baseline-adjusted time series for the parameter. If not provided, global default value is used.                          |
| timeseries_features_to_calculate | chr   |               | Comma-delimited string of features to calculate for this parameter. Overrides the features listed in default_timeseries_features_to_calculate.     |
| use_only_custom_timeseries | boolean   |  Y             | If set to TRUE, only use the custom timeseries defined in custom_timeseries - even if autogenerate_timeseries is set to TRUE.     |

Example for two parameters, one vital sign and one local lab. For any
time series generated for the vital sign, at least 50 eligible subjects
are required. For the local lab, time series with one time point only
are acceptable. These parameter-specific properties override the
corresponding global parameters (see above).

| parameter_id  | parameter_category_1 | parameter_category_2 | parameter_category_3 | parameter_name | time_point_count_min | subject_count_min | max_share_missing | generate_change_from_baseline | timeseries_features_to_calculate | use_only_custom_timeseries |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| dummy_param_1 | VS                   | NA                   | NA                   | DIABP          |                      | 50                |                   |                               |                                  | TRUE |
| dummy_param_2 | LB                   | General chemistry    | Local lab            | Glucose        | 1                    |                   |                   |                               | average;sd                       | FALSE |

### data

Data frame for the data collected during the study with one row per data
point.

Data frame columns:

| Column name      | Data type | Data required | Description                                                                                                            |
|------------------|------------------|------------------|------------------|
| subject_id       | chr       | Y             | Pointer to the subject identifier in the subjects table.                                                               |
| parameter_id     | chr       | Y             | Pointer to the parameter identifier in the parameters table.                                                           |
| timepoint_1_name | chr       | Y             | Name of the top level time point (e.g. Visit 1).                                                                       |
| timepoint_2_name | chr       |               | Name of the second level time point (e.g. “30 mins after dosing”).                                                     |
| timepoint_rank   | int       | Y             | Integer for ranking time points into chronological order. Ranking is separate for each parameter.                      |
| result           | float     | Y             | Numerical result of the parameter.                                                                                     |
| baseline         | float     |               | Possible baseline value for the parameter and the subject. This is used to calculate change-from-baseline time series. |

Example of data contents:

| subject_id | parameter_id  | timepoint_1_name | timepoint_2_name | timepoint_rank | result | baseline |
|-----------|-----------|-----------|-----------|-----------|-----------|-----------|
| subj_1     | dummy_param_1 | Visit 1          | pre-dose         | 1              | 72     | 65       |
| subj_1     | dummy_param_1 | Visit 1          | post-dose        | 2              | 81     | 65       |
| subj_1     | dummy_param_1 | Visit 2          | NA               | 3              | 69     | 65       |
| subj_1     | dummy_param_1 | Visit 3          | NA               | 4              | 70     | 65       |
| subj_1     | dummy_param_2 | Visit 3          | NA               | 1              | 85     | 84       |
| subj_2     | dummy_param_2 | Visit 3          | NA               | 1              | 102    | 96       |

### custom_timeseries

The study team might be interested in specific time series - for
instance in an oncology trial one might be interested in data collected
in the first two treatment cycles.

When time series are auto-generated in a data-driven manner, it cannot
be guaranteed that a particular time series is generated. To make sure
that a specific time series is analyzed, this can be defined with the
*custom_timeseries* data frame.

Please note: if you do not have any custom time series to process, you
still have to provide the main function a blank data frame!

Data frame columns:

| Column name     | Data type | Data required | Description                                                                                                                                          |
|------------------|------------------|------------------|------------------|
| timeseries_id   | chr       | Y             | Identifier for the the custom time series.                                                                                                           |
| parameter_id    | chr       | Y             | Pointer to the parameter identifier in the parameters table.                                                                                         |
| timepoint_combo | chr       | Y             | Semicolon-delimited string of timepoint ranks that the time series consists of. These must match values in the timepoint_rank column of the data df. |

In the example below, a custom time series has been defined for
dummy_param_1. Timepoints are defined with their respective ranks in the
*data* df.

| timeseries_id | parameter_id  | timepoint_combo |
|---------------|---------------|-----------------|
| custom_ts_1   | dummy_param_1 | 1;2;3;4         |

### custom_reference_groups

For some parameters and features, it makes no sense to compare a site to all other sites globally.
These include parameters whose values depend heavily on the underlying population
such as height and weight. With this table, the user can define parameters and
features which should be compared to other sites in the same country or the same region of the world
(e.g. the average subject heights from a site in Japan should only be compared the average subject
heights in other sites in Eastern Asia).

If a parameter-feature combination is not listed here, the comparison is done by default
to all sites globally.

| Column name     | Data type | Data required | Description                                                                                                                                          |
|------------------|------------------|------------------|------------------|
| parameter_id | chr | Y | Reference to a parameter in the parameters table |
| feature | chr | Y | Time series feature |
| ref_group | chr | Y | Scope of comparison, accepted values "country", "region", "global" |

## Output data

Output is a list of four data frames. Their details are explained in the
following.

### timeseries

Time series processed. Contains both the auto-generated and/or custom
time series.

Data frame columns:

| Column name              | Data type | Description                                                                           |
|------------------------|------------------------|------------------------|
| timeseries_id            | chr       | Unique identifier for the time series.                                                |
| parameter_id             | chr       | The parameter the time series is for.                                                 |
| baseline                 | chr       | Whether the time series values are original measurements or baseline-adjusted.        |
| timepoint_combo          | chr       | Ranks of the time points that constitute the time series. Semicolon delimited string. |
| timepoint_combo_readable | chr       | Names of the time points that constitute the time series. Semicolon-delimited string. |
| timepoint_count          | int       | Number of timepoints in the time series.                                              |

### timeseries_features

Data frame with features calculated for time series.

Data frame columns:

| Column name   | Data type | Description                                                           |
|------------------------|------------------------|------------------------|
| timeseries_id | chr       | Unique identifier of the time series. Points to the timeseries table. |
| subject_id    | chr       | Subject ID pointing to the input df subjects.                         |
| feature       | chr       | Feature code.                                                         |
| feature_value | float     | Feature value.                                                        |

### PCA_coordinates

Data frame with the top two principal components for a time series.
These are useful for the visualization of the similarity of subjects.

Data frame columns:

| Column name   | Data type | Description                                                           |
|------------------------|------------------------|------------------------|
| timeseries_id | chr       | Unique identifier of the time series. Points to the timeseries table. |
| subject_id    | chr       | Subject ID pointing to the input df subjects.                         |
| pc1           | float     | The first principal component.                                        |
| pc2           | float     | The second principal component.                                       |

### site_scores

Biasness scores for sites. The data frame contains a row for each
site/time series/feature combination.

Data frame columns:

| Column name               | Data type | Description                                                                                                                  |
|------------------------|------------------------|------------------------|
| timeseries_id             | chr       | Unique identifier of the time series. Points to the timeseries table.                                                        |
| site                      | chr       | Study site for which the score has been calculated.                                                                          |
| country                   | chr       | Site’s country.                                                                                                              |
| region                   | chr       | Region of the world the site and its country belong to.                                                                                                              |
| feature                   | chr       | Time series feature the score is for.                                                                                        |
| pvalue_kstest_logp        | float     | Negative logarithm of the raw p-value.                                                                                       |
| kstest_statistic          | float     | Test statistic of the Kolmogorov-Smirnov test.                                                                               |
| fdr_corrected_pvalue_logp | float     | Site score: negative logarithm of the multiple testing corrected p-value. FDR (False Discovery Rate) is used for correction. |
| ref_group           | chr       | Reference group used for calculating the score.                        |
| subject_count             | int       | Number of site subjects who are eligible for the time series.                                                                |
