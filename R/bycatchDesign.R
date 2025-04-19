
#---------------------------------
#Design-based estimators (new)
#----------------------------------

#Roxygen header
#'Bycatch estimation using design-based estimators
#'
#' ADD DESCRIPTION
#'
#'
#' @param setupObj  An object produced by \code{bycatchSetup_new}.
#' @param designMethods Character vector of methods to use for design based estimation. Current options are Ratio and Delta (for a delta-lognormal estimator).
#' @param designVars Specify strata that must be included in design based estimates, in order across which data should be pooled
#' @param designPooling TRUE if design-based estimates should be pooled for strata with missing data
#' @param poolTypes Type of pooling for each variable in designVars, as a character vector in the same order. Options are "all", "pooledVar" and (currently for year only) "adjacent"
#' @param pooledVar Variables to pool over for any variable with pooledVar in the previous line, as a character vector in the same order as designVars. Use NA for variables with other pooling methods.  This can be used to pool (for example) months into seasons when pooling is needed.
#' @param adjacentNum Number of adjacent years to include for adjacent pooling, as a numberical vector in the same order as designVars. NA for anything other than year.
#' @param minStrataUnit The smallest sample size in the strata defined by designVars that is acceptable, in sample units (e.g. trips)
#' @param baseDir Character. A directory to save output. Defaults to current working directory.
#' @import ggplot2 dplyr utils tidyverse
#' @importFrom stats median
#' @export
#' @keywords Fitting functions
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
#' #-------------------------------------------------
#' #Step 1. Run the datasetup function and review data inputs
#'setupObj<-bycatchSetup_new( #reorder arguments
#' obsdat = obsdatExample,
#' logdat = logdatExample,
#' yearVar = "Year",
#' obsEffort = "sampled.sets",
#' logEffort = "sets",
#' logUnsampledEffort = NULL,
#' includeObsCatch  = FALSE,
#' matchColumn = NA,
#' factorVariables = c("Year","season"),
#' numericVariables = "Value",
#' EstimateIndex = TRUE,
#' EstimateBycatch = TRUE,
#' logNum = NA,
#' sampleUnit = "trips",
#' baseDir = getwd(),
#' runName = "SimulatedExample",
#' runDescription = "Example with simulated data",
#' common = "Simulated species",
#' sp = "Genus species",
#' obsCatch = "Catch",
#' catchUnit = "number",
#' catchType = "dead discard"
#')
#'
#'-------------
#' #Step 2. Design-based estimators (with pooling)
#'bycatchDesign(
#' setupObj = setupObj,
#' designMethods = c("Ratio", "Delta"),
#' designVars = c("Year","season"),
#' designPooling = TRUE,
#' poolTypes=c("adjacent","all"),
#' pooledVar=c(NA,NA),
#' adjacentNum=c(1,NA),
#' minStrataUnit = 1
#')}


bycatchDesign_new <- function(
  setupObj = setupObj,
  designMethods = "Ratio",
  designVars = "Year",
  designPooling = FALSE,
  poolTypes = NULL,
  pooledVar = NULL,
  adjacentNum = NULL,
  minStrataUnit = 1,
  baseDir = getwd()
){






}



