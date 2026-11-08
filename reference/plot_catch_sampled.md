# Make a plot of the observed effort and catches

Make a plot of the observed effort and catches

## Usage

``` r
plot_catch_sampled(
  obsdat,
  bin_size = 1,
  effort_col = "Effort",
  catch_col = "Catch"
)
```

## Arguments

- obsdat:

  observer data

- bin_size:

  spatial bin size

- effort_col:

  column with effort data

- catch_col:

  column with catch data

## Value

Makes a plot of sampled catch and effort, with effort as tile fill and
catch as circle size
