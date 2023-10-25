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
#> # A tibble: 7,512 × 7
#>    subject_id timepoint_rank timepoint_1_name result parameter…¹ timep…² basel…³
#>    <chr>               <int> <chr[1d]>         <dbl> <chr>       <chr>   <lgl>  
#>  1 1                       1 AB                 37.6 param1      timepo… NA     
#>  2 1                       2 AC                 25.3 param1      timepo… NA     
#>  3 1                       3 AD                 NA   param1      timepo… NA     
#>  4 1                       4 AE                 32.6 param1      timepo… NA     
#>  5 1                       5 AF                 26.9 param1      timepo… NA     
#>  6 1                       6 AG                 31.2 param1      timepo… NA     
#>  7 1                       7 AH                 35.9 param1      timepo… NA     
#>  8 1                       8 AI                 35.8 param1      timepo… NA     
#>  9 1                       9 AJ                 30.8 param1      timepo… NA     
#> 10 1                      10 AK                 32.7 param1      timepo… NA     
#> # … with 7,502 more rows, and abbreviated variable names ¹​parameter_id,
#> #   ²​timepoint_2_name, ³​baseline
#> 
#> $parameters
#> # A tibble: 2 × 9
#>   parameter_id paramet…¹ param…² param…³ param…⁴ time_…⁵ subje…⁶ max_s…⁷ gener…⁸
#>   <chr>        <chr>     <chr>   <chr>   <chr>   <lgl>   <lgl>   <lgl>   <lgl>  
#> 1 param1       param1    catego… catego… catego… NA      NA      NA      NA     
#> 2 param2       param2    catego… catego… catego… NA      NA      NA      NA     
#> # … with abbreviated variable names ¹​parameter_name, ²​parameter_category_1,
#> #   ³​parameter_category_2, ⁴​parameter_category_3, ⁵​time_point_count_min,
#> #   ⁶​subject_count_min, ⁷​max_share_missing, ⁸​generate_change_from_baseline
#> 
#> $subjects
#> # A tibble: 189 × 3
#>    subject_id site  country
#>    <chr>      <chr> <chr>  
#>  1 1          1     C      
#>  2 2          1     C      
#>  3 3          1     C      
#>  4 4          1     C      
#>  5 5          1     C      
#>  6 6          1     C      
#>  7 7          1     C      
#>  8 8          2     C      
#>  9 9          2     C      
#> 10 10         2     C      
#> # … with 179 more rows
#> 
#> $custom_timeseries
#> # A tibble: 0 × 0

ls_tsoa <- process_a_study(
  data = tsoa_data$data,
  subjects = tsoa_data$subjects,
  parameters = tsoa_data$parameters,
  custom_timeseries = tsoa_data$custom_timeseries,
  timeseries_features_to_calculate = c(
    "autocorr",
    "average",
    "own_site_simil_score",
    "sd",
    "unique_value_count_relative",
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
#> # A tibble: 12 × 6
#>    timeseries_id          parameter_id baseline timepoint_combo  timep…¹ timep…²
#>    <chr>                  <chr>        <chr>    <chr>            <chr>     <dbl>
#>  1 ts_1_autogen_original  param1       original 1;2;3;4;5;6;7;8… AB_tim…      35
#>  2 ts_2_autogen_original  param1       original 1;2;3;4;5;6;7;8… AB_tim…      34
#>  3 ts_3_autogen_original  param1       original 1;2;3;4;5;6;7;8… AB_tim…      32
#>  4 ts_4_autogen_original  param1       original 1;2;3;4;5;6;7;8… AB_tim…      30
#>  5 ts_5_autogen_original  param1       original 1;2;3;4;5;6;7;8… AB_tim…      28
#>  6 ts_6_autogen_original  param1       original 1;2;3;4;5;6;7;8… AB_tim…      26
#>  7 ts_7_autogen_original  param1       original 1;2;3;4;5;6;7;8… AB_tim…      22
#>  8 ts_8_autogen_original  param2       original 1;2;3;4;5;6;7;8… AB_tim…      34
#>  9 ts_9_autogen_original  param2       original 1;2;3;4;5;6;7;8… AB_tim…      32
#> 10 ts_10_autogen_original param2       original 1;2;3;4;5;6;7;8… AB_tim…      30
#> 11 ts_11_autogen_original param2       original 1;2;3;4;5;6;7;8… AB_tim…      26
#> 12 ts_12_autogen_original param2       original 1;2;3;4;5;6;7;8… AB_tim…      22
#> # … with abbreviated variable names ¹​timepoint_combo_readable, ²​timepoint_count
#> 
#> $timeseries_features
#> # A tibble: 4,219 × 6
#>    timeseries_id         subject_id feature                featu…¹ site  country
#>    <chr>                 <chr>      <chr>                    <dbl> <chr> <chr>  
#>  1 ts_1_autogen_original 1          sd                       4.21  1     C      
#>  2 ts_1_autogen_original 1          average                 33.4   1     C      
#>  3 ts_1_autogen_original 1          unique_value_count_re…   1     1     C      
#>  4 ts_1_autogen_original 8          sd                       6.77  2     C      
#>  5 ts_1_autogen_original 8          average                 33.5   2     C      
#>  6 ts_1_autogen_original 8          unique_value_count_re…   1     2     C      
#>  7 ts_1_autogen_original 28         sd                       4.94  6     B      
#>  8 ts_1_autogen_original 28         average                 29.5   6     B      
#>  9 ts_1_autogen_original 28         unique_value_count_re…   1     6     B      
#> 10 ts_1_autogen_original 28         own_site_simil_score     0.414 6     B      
#> # … with 4,209 more rows, and abbreviated variable name ¹​feature_value
#> 
#> $PCA_coordinates
#> # A tibble: 1,076 × 4
#>    timeseries_id         subject_id    pc1     pc2
#>    <chr>                 <chr>       <dbl>   <dbl>
#>  1 ts_1_autogen_original 1          -10.9   -0.458
#>  2 ts_1_autogen_original 8           -8.61   5.69 
#>  3 ts_1_autogen_original 28          12.8   -1.92 
#>  4 ts_1_autogen_original 30         -18.1    2.89 
#>  5 ts_1_autogen_original 34          18.7   32.3  
#>  6 ts_1_autogen_original 35          51.0  -11.9  
#>  7 ts_1_autogen_original 40         -11.4  -14.3  
#>  8 ts_1_autogen_original 43          10.5   30.2  
#>  9 ts_1_autogen_original 62          17.6   -4.69 
#> 10 ts_1_autogen_original 71           4.67  -2.42 
#> # … with 1,066 more rows
#> 
#> $site_scores
#> # A tibble: 1,363 × 9
#>    timeseries_id   site  country feature pvalu…¹ kstes…² fdr_c…³ refer…⁴ subje…⁵
#>    <chr>           <chr> <chr>   <chr>     <dbl>   <dbl>   <dbl> <chr>     <int>
#>  1 ts_1_autogen_o… 1     C       sd      0.111     0.633  0      all           1
#>  2 ts_1_autogen_o… 2     C       sd      0.149     0.667  0      all           1
#>  3 ts_1_autogen_o… 6     B       sd      0.196     0.483  0      all           2
#>  4 ts_1_autogen_o… 7     C       sd      0.810     0.759  0.0452 all           2
#>  5 ts_1_autogen_o… 8     D       sd      0.287     0.767  0      all           1
#>  6 ts_1_autogen_o… 9     A       sd      0.889     0.967  0.0756 all           1
#>  7 ts_1_autogen_o… 13    B       sd      0.287     0.767  0      all           1
#>  8 ts_1_autogen_o… 15    B       sd      0.00466   0.238  0      all           3
#>  9 ts_1_autogen_o… 17    C       sd      0.412     0.833  0.0174 all           1
#> 10 ts_1_autogen_o… 18    D       sd      0.407     0.586  0.0174 all           2
#> # … with 1,353 more rows, and abbreviated variable names ¹​pvalue_kstest_logp,
#> #   ²​kstest_statistic, ³​fdr_corrected_pvalue_logp, ⁴​reference_group,
#> #   ⁵​subject_count
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
| range                       | The range of values in a time series                             |

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
