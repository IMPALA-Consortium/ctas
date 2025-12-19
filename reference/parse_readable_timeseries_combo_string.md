# parse_readable_timeseries_combo_string

Parse timepoint codes into a more human-readable format.

## Usage

``` r
parse_readable_timeseries_combo_string(
  timeseries_combo,
  this_parameter_id,
  this_timepoint_rank_to_name_mapping
)
```

## Arguments

- timeseries_combo:

  Semicolon-delimited string with time point codes.

- this_parameter_id:

  Parameter the time point combo is for.

- this_timepoint_rank_to_name_mapping:

  Mapping table from time point code to human-readable name.

## Value

Semicolon-delimited string of human-readable timepoint names.

## Author

Pekka Tiikkainen, <pekka.tiikkainen@bayer.com>
