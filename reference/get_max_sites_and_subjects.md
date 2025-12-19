# Get timpepoint that has the most sites and patients

Get timpepoint that has the most sites and patients

## Usage

``` r
get_max_sites_and_subjects(
  timepoint_ranks,
  this_time_point_count_min,
  this_max_share_missing,
  dataset,
  this_subject_count_min,
  subjects
)
```

## Arguments

- timepoint_ranks:

  a vector with all time point ranks

- this_time_point_count_min:

  Minimum number of time points a time series must have.

- this_max_share_missing:

  Maximum share of missing measuremnents a subject can have.

- dataset:

  Description of parameter x.

- this_subject_count_min:

  Minimum number of subjects per time series.

- subjects:

  Data frame with one row per study subject.

## Value

integer with optimized timepoint

## Details

is called by pick_timepoint combos. uses the same for loop to iterate
over possible timeseries and returns eligable patients and sites.
returns highest timepoint rank that has the highest number of sites with
the highest number of patients. As NA values can occur randomly patient
elegibility can vary depending on the timepoint range. It is not
necessarily the shortest range that has the most patients. This is why
we iteratively need to determine the optimal conditions using the for
loop.
