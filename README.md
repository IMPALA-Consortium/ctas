<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R-CMD-check](https://github.com/IMPALA-Consortium/ctas/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/IMPALA-Consortium/ctas/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

# Clinical Timeseries Anomaly Spotter (CTAS)

Pekka Tiikkainen, Bayer Pharmaceuticals

## Introduction

CTAS is an R package for identifying outliers and anomalies in clinical
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

library(ctas)
#> Loading required package: dplyr
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

data("ctas_data", package = "ctas")

ctas_data
#> $data
#> # A tibble: 6,684 × 7
#>    subject_id timepoint_rank timepoint_1_name result parameter…¹ timep…² basel…³
#>    <chr>               <int> <chr[1d]>         <dbl> <chr>       <chr>   <lgl>  
#>  1 1                       1 AB                 NA   param1      timepo… NA     
#>  2 1                       2 AC                 31.8 param1      timepo… NA     
#>  3 1                       3 AD                 36.4 param1      timepo… NA     
#>  4 1                       4 AE                 33.5 param1      timepo… NA     
#>  5 1                       5 AF                 34.1 param1      timepo… NA     
#>  6 1                       6 AG                 NA   param1      timepo… NA     
#>  7 1                       7 AH                 31.3 param1      timepo… NA     
#>  8 1                       8 AI                 33.2 param1      timepo… NA     
#>  9 1                       9 AJ                 35.0 param1      timepo… NA     
#> 10 1                      10 AK                 NA   param1      timepo… NA     
#> # … with 6,674 more rows, and abbreviated variable names ¹​parameter_id,
#> #   ²​timepoint_2_name, ³​baseline
#> 
#> $parameters
#> # A tibble: 2 × 11
#>   parameter_id paramet…¹ param…² param…³ param…⁴ time_…⁵ subje…⁶ max_s…⁷ gener…⁸
#>   <chr>        <chr>     <chr>   <chr>   <chr>   <lgl>   <lgl>   <lgl>   <lgl>  
#> 1 param1       param1    catego… catego… catego… NA      NA      NA      NA     
#> 2 param2       param2    catego… catego… catego… NA      NA      NA      NA     
#> # … with 2 more variables: timeseries_features_to_calculate <lgl>,
#> #   use_only_custom_timeseries <lgl>, and abbreviated variable names
#> #   ¹​parameter_name, ²​parameter_category_1, ³​parameter_category_2,
#> #   ⁴​parameter_category_3, ⁵​time_point_count_min, ⁶​subject_count_min,
#> #   ⁷​max_share_missing, ⁸​generate_change_from_baseline
#> 
#> $subjects
#> # A tibble: 162 × 4
#>    subject_id site  country region
#>    <chr>      <chr> <chr>   <chr> 
#>  1 1          AAA   AA      A     
#>  2 2          AAA   AA      A     
#>  3 3          AAA   AA      A     
#>  4 4          AAA   AA      A     
#>  5 5          AAA   AA      A     
#>  6 6          AAA   AA      A     
#>  7 7          AAB   AA      A     
#>  8 8          AAB   AA      A     
#>  9 9          AAB   AA      A     
#> 10 10         AAB   AA      A     
#> # … with 152 more rows
#> 
#> $custom_timeseries
#> # A tibble: 0 × 3
#> # … with 3 variables: timeseries_id <chr>, parameter_id <chr>,
#> #   timepoint_combo <chr>
#> 
#> $custom_reference_groups
#> # A tibble: 0 × 3
#> # … with 3 variables: parameter_id <chr>, feature <chr>, ref_group <chr>

feats <- c(
    "autocorr",
    "average",
    "own_site_simil_score",
    "sd",
    "unique_value_count_relative",
    "lof",
    "range"
  ) %>%
  paste(collapse = ";")

ls_ctas <- process_a_study(
  data = ctas_data$data,
  subjects = ctas_data$subjects,
  parameters = ctas_data$parameters,
  custom_timeseries = ctas_data$custom_timeseries,
  custom_reference_groups = ctas_data$custom_reference_groups,
  default_timeseries_features_to_calculate = feats,
  default_minimum_timepoints_per_series = 3,
  default_minimum_subjects_per_series = 3,
  default_max_share_missing_timepoints_per_series = 0.5,
  default_generate_change_from_baseline = FALSE,
  autogenerate_timeseries = TRUE
)

ls_ctas
#> $timeseries
#> # A tibble: 8 × 6
#>   timeseries_id         parameter_id baseline timepoint_combo    timep…¹ timep…²
#>   <chr>                 <chr>        <chr>    <chr>              <chr>     <dbl>
#> 1 ts_1_autogen_original param1       original 1;2;3;4;5;6;7;8;9… AB_tim…      31
#> 2 ts_2_autogen_original param1       original 1;2;3;4;5;6;7;8;9… AB_tim…      28
#> 3 ts_3_autogen_original param1       original 1;2;3;4;5;6;7;8;9… AB_tim…      24
#> 4 ts_4_autogen_original param1       original 1;2;3;4;5;6;7;8;9… AB_tim…      20
#> 5 ts_5_autogen_original param2       original 1;2;3;4;5;6;7;8;9… AB_tim…      31
#> 6 ts_6_autogen_original param2       original 1;2;3;4;5;6;7;8;9… AB_tim…      30
#> 7 ts_7_autogen_original param2       original 1;2;3;4;5;6;7;8;9… AB_tim…      26
#> 8 ts_8_autogen_original param2       original 1;2;3;4;5;6;7;8;9… AB_tim…      22
#> # … with abbreviated variable names ¹​timepoint_combo_readable, ²​timepoint_count
#> 
#> $timeseries_features
#> # A tibble: 5,706 × 7
#>    timeseries_id         subject_id feature        featur…¹ site  country region
#>    <chr>                 <chr>      <chr>             <dbl> <chr> <chr>   <chr> 
#>  1 ts_1_autogen_original 1          range           6.44    AAA   AA      A     
#>  2 ts_1_autogen_original 1          sd              1.66    AAA   AA      A     
#>  3 ts_1_autogen_original 1          unique_value_…  1       AAA   AA      A     
#>  4 ts_1_autogen_original 1          autocorr       -0.393   AAA   AA      A     
#>  5 ts_1_autogen_original 1          average        33.9     AAA   AA      A     
#>  6 ts_1_autogen_original 1          lof             0.996   AAA   AA      A     
#>  7 ts_1_autogen_original 11         range           6.84    AAC   AA      A     
#>  8 ts_1_autogen_original 11         sd              1.90    AAC   AA      A     
#>  9 ts_1_autogen_original 11         unique_value_…  1       AAC   AA      A     
#> 10 ts_1_autogen_original 11         autocorr       -0.00259 AAC   AA      A     
#> # … with 5,696 more rows, and abbreviated variable name ¹​feature_value
#> 
#> $PCA_coordinates
#> # A tibble: 822 × 4
#>    timeseries_id         subject_id    pc1     pc2
#>    <chr>                 <chr>       <dbl>   <dbl>
#>  1 ts_1_autogen_original 1          -16.6   -1.57 
#>  2 ts_1_autogen_original 11          16.1    0.857
#>  3 ts_1_autogen_original 15          27.9   -1.62 
#>  4 ts_1_autogen_original 16         -34.4    6.17 
#>  5 ts_1_autogen_original 17           7.19  -2.61 
#>  6 ts_1_autogen_original 20           6.33  -0.144
#>  7 ts_1_autogen_original 21           5.30  14.9  
#>  8 ts_1_autogen_original 26         -32.1   -3.02 
#>  9 ts_1_autogen_original 28          -1.96 -21.4  
#> 10 ts_1_autogen_original 31           1.33   5.71 
#> # … with 812 more rows
#> 
#> $site_scores
#> # A tibble: 1,765 × 10
#>    timese…¹ site  country region feature pvalu…² kstes…³ fdr_c…⁴ ref_g…⁵ subje…⁶
#>    <chr>    <chr> <chr>   <chr>  <chr>     <dbl>   <dbl>   <dbl> <chr>     <int>
#>  1 ts_1_au… AAA   AA      A      range    0.680    0.909 0.0305  global        1
#>  2 ts_1_au… AAC   AA      A      range    0.571    0.879 0.0305  global        1
#>  3 ts_1_au… AAD   AA      A      range    0.191    0.391 0.00685 global        3
#>  4 ts_1_au… ABA   AB      A      range    1.08     0.603 0.0905  global        4
#>  5 ts_1_au… ABB   AB      A      range    0.0966   0.415 0.00312 global        2
#>  6 ts_1_au… ABC   AB      A      range    0.248    0.492 0.0108  global        2
#>  7 ts_1_au… ACA   AC      A      range    0.391    0.421 0.0108  global        4
#>  8 ts_1_au… ADA   AD      A      range    0.571    0.879 0.0305  global        1
#>  9 ts_1_au… ADB   AD      A      range    0.934    0.571 0.0402  global        4
#> 10 ts_1_au… ADC   AD      A      range    0.910    0.769 0.0402  global        2
#> # … with 1,755 more rows, and abbreviated variable names ¹​timeseries_id,
#> #   ²​pvalue_kstest_logp, ³​kstest_statistic, ⁴​fdr_corrected_pvalue_logp,
#> #   ⁵​ref_group, ⁶​subject_count
```

## Input parameters

The function process_a_study has a number of parameters that need to be
provided. Some of these are data frames while others are simple
variables.

The following gives more details on each parameter. Examples are also
provided where needed.

### timeseries_features_to_calculate

Vector of features to calculate for each time series. The vector must
contain at least one of the following feature codes:

| Feature code                | Description                                                      |
|------------------------------------|------------------------------------|
| autocorr                    | Auto-correlation                                                 |
| average                     | Average value                                                    |
| own_site_simil_score        | Measure of co-clustering of time series from the same study site |
| sd                          | Standard deviation                                               |
| unique_value_count_relative | Number of unique values divided by number of values available    |

### default_minimum_timepoints_per_series

Minimum number of time points for auto-generated time series. This value
will be used for a parameter if they do not have the
time_point_count_min column defined in the parameters df.

### default_minimum_subjects_per_series

Minimum number of eligible subjects for auto-generated time series. This
value will be used for a parameter if they do not have the
subject_count_min column defined in the parameters df. The parameter
has a minimum value of two.

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
| subject_id  | chr       | Y             | Unique identifier for each subject     |
| country     | chr       | Y             | Country where the subject is enrolled. |
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

Example for two parameters, one vital sign and one local lab. For any
time series generated for the vital sign, at least 50 eligible subjects
are required. For the local lab, time series with one time point only
are acceptable. These parameter-specific properties override the
corresponding global parameters (see above).

| parameter_id  | parameter_category_1 | parameter_category_2 | parameter_category_3 | parameter_name | time_point_count_min | subject_count_min | max_share_missing | generate_change_from_baseline |
|--------|--------|--------|--------|--------|--------|--------|--------|--------|
| dummy_param_1 | VS                   | NA                   | NA                   | DIABP          |                      | 50                |                   |                               |
| dummy_param_2 | LB                   | General chemistry    | Local lab            | Glucose        | 1                    |                   |                   |                               |

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
| feature                   | chr       | Time series feature the score is for.                                                                                        |
| pvalue_kstest_logp        | float     | Negative logarithm of the raw p-value.                                                                                       |
| kstest_statistic          | float     | Test statistic of the Kolmogorov-Smirnov test.                                                                               |
| fdr_corrected_pvalue_logp | float     | Site score: negative logarithm of the multiple testing corrected p-value. FDR (False Discovery Rate) is used for correction. |
| reference_group           | chr       | Reference group for the site. “All” means that the site has been compared to all sites in the study.                         |
| subject_count             | int       | Number of site subjects who are eligible for the time series.                                                                |
