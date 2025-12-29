# Function to make design based estimates of bycatch from the ratio estimator of Pennington Delta estimator, pooling as needed for strata missing data. stratification defined by designVars, then aggregated to strataVars

Function to make design based estimates of bycatch from the ratio
estimator of Pennington Delta estimator, pooling as needed for strata
missing data. stratification defined by designVars, then aggregated to
strataVars

## Usage

``` r
getDesignEstimates(
  obsdatval,
  logdatval,
  strataVars,
  designVars = NULL,
  designPooling,
  minStrataUnit = 1,
  startYear,
  poolingSum = NULL,
  includePool = NULL
)
```

## Arguments

- obsdatval:

  Value

- logdatval:

  Value

- strataVars:

  Value

- designVars:

  Value

- designPooling:

  Value

- minStrataUnit:

  Value

- startYear:

  Value

- poolingSum:

  Value

- includePool:

  Value
