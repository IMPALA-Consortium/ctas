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
#> # A tibble: 6,890 × 7
#>    subject_id timepoint_rank timepoint_1_name result parameter_id
#>    <chr>               <int> <chr[1d]>         <dbl> <chr>       
#>  1 1                       1 AB                 38.1 param1      
#>  2 1                       2 AC                 36.9 param1      
#>  3 1                       3 AD                 38.5 param1      
#>  4 1                       4 AE                 39.8 param1      
#>  5 1                       5 AF                 40.6 param1      
#>  6 1                       6 AG                 36.7 param1      
#>  7 1                       7 AH                 38.2 param1      
#>  8 1                       8 AI                 40.0 param1      
#>  9 1                       9 AJ                 37.0 param1      
#> 10 1                      10 AK                 35.1 param1      
#> # ℹ 6,880 more rows
#> # ℹ 2 more variables: timepoint_2_name <chr>, baseline <lgl>
#> 
#> $parameters
#> # A tibble: 2 × 11
#>   parameter_id parameter_name parameter_category_1 parameter_category_2
#>   <chr>        <chr>          <chr>                <chr>               
#> 1 param1       param1         category 1           category 2          
#> 2 param2       param2         category 1           category 2          
#> # ℹ 7 more variables: parameter_category_3 <chr>, time_point_count_min <lgl>,
#> #   subject_count_min <lgl>, max_share_missing <lgl>,
#> #   generate_change_from_baseline <lgl>,
#> #   timeseries_features_to_calculate <lgl>, use_only_custom_timeseries <lgl>
#> 
#> $subjects
#> # A tibble: 177 × 4
#>    subject_id site  country region
#>    <chr>      <chr> <chr>   <chr> 
#>  1 1          AAA   AA      A     
#>  2 2          AAA   AA      A     
#>  3 3          AAA   AA      A     
#>  4 4          AAB   AA      A     
#>  5 5          AAB   AA      A     
#>  6 6          AAB   AA      A     
#>  7 7          AAC   AA      A     
#>  8 8          AAC   AA      A     
#>  9 9          AAC   AA      A     
#> 10 10         AAC   AA      A     
#> # ℹ 167 more rows
#> 
#> $custom_timeseries
#> # A tibble: 0 × 3
#> # ℹ 3 variables: timeseries_id <chr>, parameter_id <chr>, timepoint_combo <chr>
#> 
#> $custom_reference_groups
#> # A tibble: 0 × 3
#> # ℹ 3 variables: parameter_id <chr>, feature <chr>, ref_group <chr>

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
#> # A tibble: 10 × 6
#>    timeseries_id    parameter_id baseline timepoint_combo timepoint_combo_read…¹
#>    <chr>            <chr>        <chr>    <chr>           <chr>                 
#>  1 ts_1_autogen_or… param1       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  2 ts_2_autogen_or… param1       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  3 ts_3_autogen_or… param1       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  4 ts_4_autogen_or… param1       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  5 ts_5_autogen_or… param1       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  6 ts_6_autogen_or… param2       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  7 ts_7_autogen_or… param2       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  8 ts_8_autogen_or… param2       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#>  9 ts_9_autogen_or… param2       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#> 10 ts_10_autogen_o… param2       original 1;2;3;4;5;6;7;… AB_timepoint name 2;A…
#> # ℹ abbreviated name: ¹​timepoint_combo_readable
#> # ℹ 1 more variable: timepoint_count <dbl>
#> 
#> $timeseries_features
#> # A tibble: 7,902 × 7
#>    timeseries_id         subject_id feature   feature_value site  country region
#>    <chr>                 <chr>      <chr>             <dbl> <chr> <chr>   <chr> 
#>  1 ts_1_autogen_original 1          range            5.45   AAA   AA      A     
#>  2 ts_1_autogen_original 1          sd               1.40   AAA   AA      A     
#>  3 ts_1_autogen_original 1          unique_v…        1      AAA   AA      A     
#>  4 ts_1_autogen_original 1          autocorr         0.120  AAA   AA      A     
#>  5 ts_1_autogen_original 1          average         38.2    AAA   AA      A     
#>  6 ts_1_autogen_original 1          lof              1.06   AAA   AA      A     
#>  7 ts_1_autogen_original 4          range           14.6    AAB   AA      A     
#>  8 ts_1_autogen_original 4          sd               3.70   AAB   AA      A     
#>  9 ts_1_autogen_original 4          unique_v…        1      AAB   AA      A     
#> 10 ts_1_autogen_original 4          autocorr        -0.0317 AAB   AA      A     
#> # ℹ 7,892 more rows
#> 
#> $PCA_coordinates
#> # A tibble: 1,137 × 4
#>    timeseries_id         subject_id      pc1    pc2
#>    <chr>                 <chr>         <dbl>  <dbl>
#>  1 ts_1_autogen_original 1          -37.7      1.07
#>  2 ts_1_autogen_original 4           14.6      2.67
#>  3 ts_1_autogen_original 10          47.6     -2.96
#>  4 ts_1_autogen_original 12          35.0      1.42
#>  5 ts_1_autogen_original 14           0.0945 -34.7 
#>  6 ts_1_autogen_original 16           2.97    -9.54
#>  7 ts_1_autogen_original 17          -8.90     3.78
#>  8 ts_1_autogen_original 19         -29.5      2.09
#>  9 ts_1_autogen_original 22         -12.3      1.59
#> 10 ts_1_autogen_original 24           6.02    -5.75
#> # ℹ 1,127 more rows
#> 
#> $site_scores
#> # A tibble: 2,221 × 10
#>    timeseries_id         site  country region feature pvalue_kstest_logp
#>    <chr>                 <chr> <chr>   <chr>  <chr>                <dbl>
#>  1 ts_1_autogen_original AAA   AA      A      range               0.490 
#>  2 ts_1_autogen_original AAB   AA      A      range               0.0130
#>  3 ts_1_autogen_original AAC   AA      A      range               0.418 
#>  4 ts_1_autogen_original AAD   AA      A      range               1.24  
#>  5 ts_1_autogen_original AAE   AA      A      range               0.374 
#>  6 ts_1_autogen_original AAF   AA      A      range               0.0324
#>  7 ts_1_autogen_original AAG   AA      A      range               0.977 
#>  8 ts_1_autogen_original ABA   AB      A      range               1.17  
#>  9 ts_1_autogen_original BAA   BA      B      range               0.283 
#> 10 ts_1_autogen_original BAB   BA      B      range               1.77  
#> # ℹ 2,211 more rows
#> # ℹ 4 more variables: kstest_statistic <dbl>, fdr_corrected_pvalue_logp <dbl>,
#> #   ref_group <chr>, subject_count <int>
```

## Input parameters

The function process_a_study has a number of parameters that need to be
provided. Some of these are data frames while others are simple
variables.

The following gives more details on each parameter. Examples are also
provided where needed.

### timeseries_features_to_calculate

Semicolon-delimited string of time series features to calculate. The
list must contain at least one of the following feature codes:

| Feature code                | Description                                                                                                          |
|------------------------------------|------------------------------------|
| autocorr                    | Auto-correlation                                                                                                     |
| average                     | Average value                                                                                                        |
| own_site_simil_score        | Measure of co-clustering of time series from the same study site                                                     |
| sd                          | Standard deviation                                                                                                   |
| unique_value_count_relative | Number of unique values divided by number of values available                                                        |
| range                       | Maximum difference between two time points in a series                                                               |
| lof                         | Local Outlier Factor. Values around one are inliers while the larger the value, more of an outlier the timeseries is |

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

| Column name | Data type | Data required | Description                                                                                                           |
|------------------|------------------|------------------|------------------|
| subject_id  | chr       | Y             | Unique identifier for each subject                                                                                    |
| country     | chr       | Y             | Country where the subject is enrolled.                                                                                |
| site        | chr       | Y             | Name of the study site.                                                                                               |
| region      | chr       | \(Y\)         | Region of the world for the country. Required if the custom reference table has “region” defined for any of the rows. |

### parameters

Data frame with one row for each parameter that should be processed.
Parameters are, for example, laboratory assays and vital sign
measurements done during the course of a study.

Data frame columns:

| Column name                      | Data type | Data required | Description                                                                                                                                       |
|------------------|------------------|------------------|------------------|
| parameter_id                     | chr       | Y             | Parameter identifier                                                                                                                              |
| parameter_category_1             | chr       |               | Top level category of the parameter (e.g. LB for a laboratory assay).                                                                             |
| parameter_category_2             | chr       |               | Second level category of the parameter (e.g. General chemistry for a laboratory assay).                                                           |
| parameter_category_3             | chr       |               | Third level category of the parameter (e.g. Local lab for a laboratory assay).                                                                    |
| parameter_name                   | chr       | Y             | Parameter name                                                                                                                                    |
| time_point_count_min             | int       |               | Minimum number of time points required for auto-generated time series. If not provided, global default value is used.                             |
| subject_count_min                | int       |               | Minimum number of eligible subjects required for auto-generated time series. If not provided, global default value is used.                       |
| max_share_missing                | float     |               | Maximum number of missing data points an eligible subject can have for auto-generated time series. If not provided, global default value is used. |
| generate_change_from_baseline    | boolean   |               | Whether to also generate baseline-adjusted time series for the parameter. If not provided, global default value is used.                          |
| timeseries_features_to_calculate | chr       |               | A comma-delimited list of time series features to calculate for the parameter. If not provided, global default value is used.                     |
| use_only_custom_timeseries       | boolean   |               | Whether to only use the custom timeseries for the parameter (i.e. ignore any autogenerated time series).                                          |

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

### custom_reference_groups

For some parameters, it might make more sense to compare a site to other
sites in the same country or region than globally. This can be a case
for subject weight or height for example where it is better to compare,
say, the average height or weight of Japanese sites to other East Asian
sites. If compared globally, most Japanese sites would get flagged since
in general Japanese are lighter and shorter than people elsewhere.

Data frame columns:

| Column name  | Data type | Data required | Description                                                                                                                                                                             |
|------------------|------------------|------------------|------------------|
| parameter_id | chr       | Y             | Pointer to the parameter identifier in the parameters table.                                                                                                                            |
| feature      | chr       | Y             | Time series feature for which to use the custom reference group.                                                                                                                        |
| ref_group    | chr       | Y             | “country” if the parameter features should only be compared to sites from the same country, “region” if the comparison should be made to other sites from the same region of the world. |

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
