# Generate predicted total bycatch with standard errors and confidence intervals of predictions (computationally efficient version with help from GitHub copilot)

Avoids forming the n x n observation-level covariance matrix A %*% Sigma
%*% t(A) by instead computing the equivalent scalar variance as
crossprod(v, Sigma_beta %*% v) where v = t(A) %*% deriv is a p-length
vector (p = number of model parameters). This reduces memory from O(n^2)
to O(np) and computation from O(n^2 p) to O(np + p^2), making it
feasible for large n.

## Usage

``` r
makePredictionsDeltaVarFast(
  modfit1,
  modfit2 = NULL,
  newdat,
  modtype,
  obsdatval,
  includeObsCatch,
  predictionGroups = NULL,
  CIval,
  printOutput = TRUE,
  catchType,
  common,
  shortName,
  dirname,
  run,
  modelScenario,
  startYear
)
```

## Arguments

- modfit1:

  Fitted model object (binomial component for delta models, otherwise
  full model)

- modfit2:

  Fitted positive model object (Delta-Lognormal or Delta-Gamma only);
  NULL otherwise

- newdat:

  Prediction data frame (logbook data)

- modtype:

  Character string naming the model type

- obsdatval:

  Observer data frame (for includeObsCatch)

- includeObsCatch:

  Logical; whether to add observed catch to predictions

- predictionGroups:

  Character vector of stratification variable names

- CIval:

  Numeric; alpha level for confidence intervals (e.g. 0.05)

- printOutput:

  Logical; whether to write CSV outputs

- catchType:

  Value

- common:

  Value

- shortName:

  Value

- dirname:

  Value

- run:

  Value

- modelScenario:

  Value

- startYear:

  Value

## Details

For delta models (Delta-Lognormal, Delta-Gamma), the delta method is
applied to both component models jointly,with observation-level residual
variance computed from the Lo et al. (1992) variance-of-a-product
formula.
