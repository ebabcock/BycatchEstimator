
#---------------------------------
#Model-based estimators (new)
#----------------------------------

#Roxygen header
#'Bycatch estimation using model-based estimators
#'
#'Produces model-based estimates of bycatch and annual abundance index
#'
#'
#' @param setupObj  An object produced by \code{bycatchSetup_new}.
#' @param complexModel Specify as stats::formula. Specify the most complex and simplest model to be considered. The code will find compare all intermediate models using information criteria. Include only fixed effects.
#' @param simpleModel Specify as stats::formula. This model includes all variables tha must be in the final bycatch estimation model
#' @param indexModel Specify as stats::formula. Use indexModel to specify which strata to keep separate in calculating abundance indices.
#' @param modelTry  Specify which observation error models to try. Options are: "Binomial", "Normal","Lognormal", "Delta-Lognormal", and "Delta-Gamma", for models using the lm and glm functions, "NegBin" for Negative binomial using glm.nb in the MASS library, "Tweedie" for Tweedie GLM from the cpglm function in the cplm library, and "TMBbinomial","TMBnormal", "TMBlognormal", "TMBdelta-Lognormal","TMBdelta-Gamma", "TMBnbinom1", "TMBnbinom2", and "TMBtweedie" for the corresponding models from the glmmTMB library. Binomial or TMBbinomial will be run automatically as part of the delta models if any of them are selected. @param obsdat Observer data set
#' @param randomEffects Character vector. Random effects that should be included in all non-delta and binomial models, as a character vector in (e.g. "Year:area" to include Year:area as a random effect). Null if none. Note that random effects will be included in all models. The code will not evaluate whether they should be included.
#' @param randomEffects2 Character vector. Random effects that should be included in the positive catch component of delta models, as a character vector in (e.g. "Year:area" to include Year:area as a random effect). Null if none. Note that random effects will be included in all models. The code will not evaluate whether they should be included.
#' @param selectCriteria Character. Model selection criteria. Options are AICc, AIC and BIC
#' @param DoCrossValidation Specify whether to run a 10 fold cross-validation (TRUE or FALSE). This may not work with a small or unbalanced dataset
#' @param DredgeCrossValidation DredgeCrossValidation specifies whether to use information criteria to find the best model in cross validation, using the dredge function, or just keep the same model formula. Do not use dredge for very large datasets, as the run will be slow.
#' @param ResidualTest Logical. Specify whether to exclude models that fail the DHARMa residuals test.
#' @param CIval Specify confidence interval for total bycatch estimates. Should be the alpha level, e.g. 0.05 for 95%
#' @param VarCalc Character. Options are: "Simulate","DeltaMethod", or "None". Variance calculation method. Simulate will not work with a large number of sample units in the logbook data. The delta method for variance calculation is not implemented for the delta-lognormal or delta-gamma methods.
#' @param useParallel Logical. Whether to conduct the analysis using parallel processing. Only initialized when more that two cores are available.
#' @param nSims Number of simulations used to calculate confidence intervals. Ignored if \code{VarCalc} set to "None"
#' @param baseDir Character. A directory to save output. Defaults to current working directory.
#' @param plotValidation Logical. Validation. If you have true values of the total bycatch (for example in a simulation study). Make PlotValidation true and fill out the rest of the specification.
#' @param trueVals The data set that contains the true simulated total catches by year.
#' @param trueCols The column of the true simulated catches that contains true bycatch by year
#' @param doReport Logical. Create a markdown report of the analysis
#' @import MuMIn ggplot2 parallel dplyr doParallel foreach utils tidyverse parallelly
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
#' #Step 2. Model fitting
#'bycatchFit(
#' setupObj = setupObj,
#' complexModel = formula(y~(Year+season)^2),
#' simpleModel = formula(y~Year),
#' indexModel = formula(y~Year),
#' modelTry = c("Delta-Lognormal","TMBnbinom2"),

#')}
