# calculate_autocorrelation

Calculates the autocorrelation feature for a time series. Default lag is
one time point.

## Usage

``` r
calculate_autocorrelation(this_timeseries, lag = 1)
```

## Arguments

- this_timeseries:

  Measurements of a time series for an individual subject ranked by
  time.

- lag:

  Lag parameter for auto-correlation. Default lag is one.

## Value

Autocorrelation coefficient.
