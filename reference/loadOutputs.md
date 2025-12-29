# loadOutputs

Reads in all the R objects created during runs of the bycatchEstimator
for further analysis

## Usage

``` r
loadOutputs(
  baseDir = getwd(),
  runName,
  runDate = Sys.Date(),
  designScenarios = NULL,
  modelScenarios = NULL
)
```

## Arguments

- baseDir:

  The base directory for the runs, same as bycatchSetup.

- runName:

  The run name, same as bycatchSetup.

- runDate:

  The date when the model was run. Defaults to current date, but can be
  set to read in models previously run.

- designScenarios:

  Character vector of designScenario values from original run. NULL to
  read in no design-based results.

- modelScenarios:

  Character vector of modelScenario values from original run.

## Value

Returns a list with the setupObj from the specified run, a list called
designObjList which contains the design-based model inputs and outputs
for each designScenario, modelObjList, which is the same for the models
in modelScenarios, and a data frame called allYearEstimates which is the
annual estimates across all design-based and model-based scenarios in a
format suitable for ggplot.
