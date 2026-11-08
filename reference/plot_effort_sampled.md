# Make plot of effort and sampling distribution

Make plot of effort and sampling distribution

## Usage

``` r
plot_effort_sampled(logdat, obsdat, bin_size = 1, effort_col = "Effort")
```

## Arguments

- logdat:

  logbook data

- obsdat:

  observer data

- bin_size:

  bin size for spatial aggregation (in degrees)

- effort_col:

  name of effort column

## Value

Makes a plot of the spatial distribution of total effort and sampled
effort, with points sized by fraction sampled.
