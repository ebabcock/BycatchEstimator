# Function to make design based estimates of bycatch from the ratio estimator of Pennington Delta estimator, pooling as needed for strata missing data. All estimates are first calculated at the stratification defined by designVars, then aggregated to strataVars. Use outGroups to include any grouping Variables that are not in designVars.

Function to make design based estimates of bycatch from the ratio
estimator of Pennington Delta estimator, pooling as needed for strata
missing data. All estimates are first calculated at the stratification
defined by designVars, then aggregated to strataVars. Use outGroups to
include any grouping Variables that are not in designVars.

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
  includePool = NULL,
  outGroups = NULL
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

- outGroups:

  Value
