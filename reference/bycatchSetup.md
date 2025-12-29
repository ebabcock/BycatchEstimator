# Bycatch estimation data setup

Sets global conditions, makes a preliminary data summary and data checks
(make plots and tables of the data).

## Usage

``` r
bycatchSetup(
  obsdat,
  logdat,
  yearVar,
  obsEffort,
  logEffort,
  obsCatch,
  catchUnit,
  catchType,
  logNum = NA,
  sampleUnit,
  factorVariables,
  numericVariables,
  EstimateBycatch = TRUE,
  baseDir = getwd(),
  runName,
  runDescription,
  common,
  sp,
  reportType = "html"
)
```

## Arguments

- obsdat:

  Observer data set

- logdat:

  Logbook data set

- yearVar:

  Character. The column name of the year variable in `obsdat` and
  `logdat`. Both input files must contain the same variable name for
  year.

- obsEffort:

  Character. The column name of the effort variable in `obsdat`. This
  variable must have the same effort units as `logEffort`

- logEffort:

  Character. The column name of the effort variable in `logdat`.
  Optional and only used when estimating bycatch. This variable must
  have the same effort units as `obsEffort`

- obsCatch:

  Character vector. The name of the column(s) in `obsdat` that contain
  catch. If it is a vector, order of variable names must follow the same
  order as names provided in `common` and `sp`

- catchUnit:

  Character vector. Give units of catch (e.g., number) to go in plot
  labels. Must be a vector of the same length as `sp`

- catchType:

  Character vector. Give type of catch (e.g., dead discards) to go in
  plot labels. Must be a vector of the same length as `sp`

- logNum:

  Character vector. The name of the column in `logdat` that gives the
  number of sample units (e.g., trips or sets). If the logbook data is
  not aggregated (i.e. each row is a sample unit) set value to NA

- sampleUnit:

  Character. What is the sample unit in `logdat`? e.g. sets or trips.

- factorVariables:

  Character vector. Specify which variables should be interpreted as
  categorical, ensuring factor format on these variables. These
  variables must have identical names and factor levels in `obsdat` and
  `logdat`

- numericVariables:

  Character vector. Specify which variables should be interpreted as
  numeric. These variables must have identical names in `obsdat` and
  `logdat`. If there are no numeric variables, set numericVariables=NA.

- EstimateBycatch:

  Logical. Defaults to TRUE. If TRUE, you must provide logbook data or
  some other source of total effort to `logdat`. FALSE will produced
  data summaries of `obsdat` only.

- baseDir:

  Character. A directory to save output. Defaults to current working
  directory.

- runName:

  Characer. Give a name to the run, which will be used to set up a
  directory for the outputs

- runDescription:

  Character. Brief summary of the run, which will be used to set up a
  directory for the outputs

- common:

  Character vector. Provide a common name for the species used in output
  filess. Can be a vector of names to do multiple species at the same
  time.

- sp:

  Character vector. Provide a scientific name for the species used in
  output files. Can be a vector of names to do multiple species at the
  same time

- reportType:

  Character. Choose type of report to be produced. Options are html
  (default), pdf or both.

## Examples

``` r
if (FALSE) { # \dontrun{
library(BycatchEstimator)
setupObj<-bycatchSetup(
obsdat = obsdatExample,
logdat = logdatExample,
yearVar = "Year",
obsEffort = "sampled.sets",
logEffort = "sets",
obsCatch = "Catch",
catchUnit = "number",
catchType = "dead discard",
logNum = NA,
sampleUnit = "trips",
factorVariables = c("Year","season"),
numericVariables = NA,
EstimateBycatch = TRUE,
baseDir = getwd(),
runName = "SimulatedExample",
runDescription = "Example with simulated data",
common = "Simulated species",
sp = "Genus species",
reportType = "html"
)} # }
```
