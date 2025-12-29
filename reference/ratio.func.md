# Basic ratio estimator with variance (Cochran)

Output is mean and standard error of bycatch by stratum and the total
bycatch with SE. Assumes unobserved strata have zero catch

## Usage

``` r
ratio.func(x, y, g, X, N, G)
```

## Arguments

- x:

  x, y and g are vectors giving the effort/catch, bycatch and stratum of
  each observed sample unit.

- y:

  x, y and g are vectors giving the effort/catch, bycatch and stratum of
  each observed sample unit.

- g:

  x, y and g are vectors giving the effort/catch, bycatch and stratum of
  each observed sample unit.

- X:

  X is the total effort/catch by stratum

- N:

  N is the total number of sample units by stratum, if available,
  otherwise total effort

- G:

  Value
