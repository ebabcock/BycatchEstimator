
#---------------------------------
#Design-based estimators (new)
#----------------------------------

#Roxygen header
#'Bycatch estimation using design-based estimators
#'
#' Produces estimates of bycatch using design-based ratio estimator and delta estimator, with the option of pooling across stratification variables.
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

  #unpack setup obj
  obsdat<-logdat<-yearVar<-obsEffort<-logEffort<-obsCatch<-catchUnit<-catchType<-
    logNum<-sampleUnit<-factorVariables<-numericVariables<-
    logUnsampledEffort<-includeObsCatch<-matchColumn<-
    EstimateIndex<-EstimateBycatch<-
    baseDir<-runName<-runDescription<-common<-sp<-NULL

  dat<-numSp<-yearSum<-allVarNames<-startYear<-strataSum<-run<-NULL

  for(r in 1:NROW(setupObj$bycatchInputs)) assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]]) #assign values of bycatchInputs to each element
  for(r in 1:NROW(setupObj$bycatchOutputs)) assign(names(setupObj$bycatchOutputs)[r],setupObj$bycatchOutputs[[r]])

  if(designPooling & length(pooledVar[!is.na(pooledVar)]>0)) temp2<-pooledVar[!is.na(pooledVar)] else temp2<-NULL

  #Make R objects to store analysis
  poolingSum<-list()
  includePool<-list()
  yearSumGraph<-list()

  # spp loop
  for(run in 1:numSp) {
  #Make annual summary - don't think we need this because yearSum has been saved in data setup
  yearSum[[run]]<-MakeSummary(
    obsdatval = dat[[run]],
    logdatval = logdat,
    strataVars = "Year", #argument within MakeSummary function, not to be changed
    EstimateBycatch = EstimateBycatch,
    startYear = startYear
  )

  if(("Ratio" %in% designMethods | "Delta" %in% designMethods) & EstimateBycatch) { #is EstimateBycatch needed?
    if(designPooling) {
      temp<-getPooling(obsdatval= dat[[run]],
                       logdatval=logdat,
                       minStrataUnit=minStrataUnit,
                       designVars=designVars,
                       pooledVar=pooledVar,
                       poolTypes=poolTypes,
                       adjacentNum=adjacentNum)
      poolingSum[[run]]<-temp[[1]]
      write.csv(poolingSum[[run]],paste0(dirname[[run]],common[run],catchType[run],"Pooling.csv"), row.names = FALSE)
      includePool[[run]]<-temp[[2]]
    } else  {
      poolingSum[[run]]<-NULL
      includePool[[run]]<-NULL
    }

    temp<-getDesignEstimates(obsdatval = dat[[run]],
                             logdatval = logdat,
                             strataVars = "Year",
                             designVars = designVars,
                             designPooling = designPooling,
                             minStrataUnit = minStrataUnit,
                             startYear = startYear,
                             poolingSum = poolingSum[[run]],
                             includePool= includePool[[run]]
    )
    yearSum[[run]]<-left_join(yearSum[[run]],temp,by="Year")
    write.csv(temp,
              paste0(dirname[[run]],common[run],catchType[run],"DesignYear.csv"), row.names = FALSE)

    #And design based stratification
    temp<-getDesignEstimates(obsdatval = dat[[run]],
                             logdatval = logdat,
                             strataVars = designVars,
                             designVars = designVars,
                             designPooling = designPooling,
                             minStrataUnit = minStrataUnit,
                             startYear = startYear,
                             poolingSum = poolingSum[[run]],
                             includePool= includePool[[run]]
    )
    write.csv(temp,
              paste0(dirname[[run]],common[run],catchType[run],"DesignStrata.csv"), row.names = FALSE)
  }

  write.csv(yearSum[[run]],
            paste0(dirname[[run]],common[run],catchType[run],"DataSummary.csv"), row.names = FALSE) #does this contain design-based estimators as well?

  #if(EstimateBycatch) { #probably don't need this if statement
    x<-list("Unstratified ratio"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$Cat,Total.se=.data$Cse))
    if("Ratio" %in% designMethods)
      x=c(x,list("Ratio"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$ratioMean,Total.se=.data$ratioSE)))
    if("Delta" %in% designMethods)
      x=c(x,list("Design Delta"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$deltaMean,Total.se=.data$deltaSE)))

    yearSumGraph[[run]]<-bind_rows(x,.id="Source")     %>%
      mutate(TotalVar=.data$Total.se^2,Total.cv=.data$Total.se/.data$Total,
             Total.mean=NA,TotalLCI=.data$Total-1.96*.data$Total.se,TotalUCI=.data$Total+1.96*.data$Total.se) %>%
      mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI))

    #Calculations at level of simple model - isn't this being produced already in bycatchsetup?
    strataSum[[run]]<-MakeSummary(
      obsdatval = dat[[run]],
      logdatval = logdat,
      strataVars = unique(c("Year",requiredVarNames)), #replace by factorVariables?
      EstimateBycatch = EstimateBycatch,
      startYear = startYear
    )
    write.csv(strataSum[[run]],
              paste0(dirname[[run]],common[run],catchType[run],"StrataSummary.csv"), row.names = FALSE)
  #}

  } #close loop for each sp


  # list of outputs to be saved in rds file - what inputs and outputs to save?
  # any outputs from bycatchSetup to carry over here?

  #Create output list
  output<-list(

    designInputs = list(
      designMethods = designMethods,
      designVars = designVars,
      designPooling = designPooling,
      poolTypes=poolTypes,
      pooledVar=pooledVar,
      adjacentNum=adjacentNum,
      minStrataUnit = minStrataUnit,
      baseDir
    ),

    designOutputs = list(
      yearSum = yearSum,
      yearSumGraph = yearSumGraph,
      strataSum = strataSum,
      poolingSum = poolingSum,
      includePool = includePool
    )

  )


  # markdown report with results





} #close main function



