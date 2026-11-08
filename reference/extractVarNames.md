# Extract variable names from a model object, including those inside smooths and polynomials

Extract variable names from a model object, including those inside
smooths and polynomials

## Usage

``` r
extractVarNames(model1)
```

## Arguments

- model1:

  a model formula or fitted model object (e.g., from lm, glm, gam,
  glmmTMB)

## Value

A character vector of unique variable names used in the model, including
those inside smooths (s(), te(), ti(), t2()), polynomials (poly()), and
transformations (I()).
