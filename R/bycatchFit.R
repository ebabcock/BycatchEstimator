
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
#' @param reportType Character. Choose type of data checks report to be produced. Options are pdf, html or both.
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
#' randomEffects = NULL,
#' randomEffects2 = NULL,
#' selectCriteria = "BIC",
#' DoCrossValidation = TRUE,
#' DredgeCrossValidation = FALSE,
#' ResidualTest = FALSE,
#' CIval = 0.05,
#' VarCalc = "Simulate",
#' useParallel = TRUE,
#' nSims = 1000,
#' baseDir = getwd(),
#' plotValidation = FALSE,
#' trueVals = NULL,
#' trueCols = NULL,
#' doReport = TRUE
#')}
#'


bycatchFit_new<-function(
  setupObj,
  complexModel,
  simpleModel,
  indexModel = NULL,
  modelTry = c("Delta-Lognormal","Delta-Gamma","TMBnbinom1","TMBnbinom2","TMBtweedie"),
  randomEffects=NULL,
  randomEffects2=NULL,
  selectCriteria = "BIC",
  DoCrossValidation = FALSE,
  DredgeCrossValidation = FALSE,
  ResidualTest = TRUE,
  CIval = 0.05,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 10,
  baseDir = getwd(),
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL,
  doReport = TRUE
){

  # unpack setup object
#do we need all terms from data setup object? are terms in model fit that are named the same and will cause issues?
#review if EstimateIndex and EstimateBycatch should actually be here or in data setup?




  NumCores<-parallelly::availableCores()  #Check if machine has multiple cores for parallel processing
  #Make sure there are multiple cores to use Parallel processing
  if(NumCores<=1) useParallel=FALSE

  #Check that all models in modelTry are valid
  if(!all(modelTry %in% c("Tweedie","Lognormal","Delta-Lognormal","Delta-Gamma", "TMBnbinom1","TMBlognormal",
                          "TMBnbinom2","TMBtweedie","Normal","Binomial","NegBin", "TMBgamma","Gamma",
                          "TMBbinomial","TMBnormal","TMBdelta-Lognormal","TMBdelta-Gamma") ))
    stop(paste("Model requested in modelTry not available"))

  #Make sure binomial is included if either of the delta models is
  if(("Delta-Lognormal" %in% modelTry |"Delta-Gamma" %in% modelTry) & !"Binomial" %in% modelTry)
    modelTry<-c("Binomial",modelTry)
  if(("TMBdelta-Lognormal" %in% modelTry |"TMBdelta-Gamma" %in% modelTry) & !"TMBbinomial" %in% modelTry)
    modelTry<-c("TMBbinomial",modelTry)

  #If there are any random effects, all fitting will be done in glmmTMB
  if(!is.null(randomEffects) | !is.null(randomEffects2)) {
    modelTry<-case_when(modelTry=="Binomial" ~"TMBbinomial",
                        modelTry=="Normal" ~"TMBnormal",
                        modelTry=="Tweedie" ~"TMBtweedie",
                        modelTry=="Gamma"~"TMBgamma",
                        modelTry=="Delta-Lognormal"~"TMBdelta-Lognormal",
                        modelTry=="Delta-Gamma"~"TMBdelta-Gamma",
                        modelTry=="Lognormal"~"TMBlognormal",
                        modelTry=="NegBin"~"TMBnbinom2",
                        grepl("TMB",modelTry)~modelTry)
    modelTry<-unique(modelTry)
  }

  requiredVarNames<-as.vector(getAllTerms(simpleModel))
  allVarNames<-as.vector(getAllTerms(complexModel))
  allVarNames<-allVarNames[grep(":",allVarNames,invert=TRUE)] #filter out interaction terms (terms that contain colon :)
  allVarNames<-allVarNames[grep("I(*)",allVarNames,invert=TRUE)] #filter out transformed terms like polynomials (I means interpret as is)
  allVarNames<-unique(c("Year",allVarNames))  #EAB 2/18/2025

  if(!is.null(randomEffects)) temp<-unlist(strsplit(randomEffects,":")) else temp<-NULL # extract random effects terms where it finds colon
  if(!is.null(randomEffects2)) temp<-c(temp,unlist(strsplit(randomEffects2,":"))) else temp<-NULL

  if(EstimateIndex) {
    indexVarNames<-as.vector(getAllTerms(indexModel))
    if(!"Year" %in% indexVarNames) indexVarNames<-c("Year",indexVarNames)
  } else indexVarNames=NULL

  if(EstimateBycatch) {
    #Add interaction random effects to logdat for ease of prediction
    randomInteractions<-unique(c(randomEffects,randomEffects2))
    for(i in which(grepl(":",randomInteractions))) {
      varval<-strsplit(randomInteractions[i],":")[[1]]
      x<-data.frame(logdat[,varval])
      x$newvar<-base::paste(x[,1],x[,2],sep=":")
      names(x)[3]<-randomInteractions[i]
      logdat<-bind_cols(logdat,select(x,randomInteractions[i]))
    }
    #Add stratum designation and check sample size in strata
    if(length(requiredVarNames) > 1) {
      logdat$strata<-apply( logdat[ , requiredVarNames ] , 1 , paste , collapse = "-" )
    }
    if(length(requiredVarNames)==1) {
      logdat$strata <- pull(logdat,var=requiredVarNames)
    }
    if(length(requiredVarNames)==0) {
      logdat$strata <- rep(1,nrow(logdat))
    }
    if(max(tapply(logdat$SampleUnits,logdat$strata,sum)) > 100000) {
      print("Inadvisable to calculate variance of estimators for such large number of logbook sample units")
    }
  }

  #indexDat for making index
  if(EstimateIndex) {
    indexDat<-distinct_at(obsdat,vars(all_of(indexVarNames)),.keep_all=TRUE) %>%
      arrange(Year) %>%
      mutate(Effort=1)
    temp<-allVarNames[allVarNames != "Year"]
    for(i in 1:length(temp)) {
      if(!temp[i] %in% indexVarNames) {
        if(is.numeric(pull(obsdat,!!temp[i])))
          indexDat[,temp[i]]<-median(pull(obsdat,!!temp[i]),na.rm=TRUE) else
            indexDat[,temp[i]]<-mostfreqfunc(obsdat[,temp[i]])
      }
    }
  } else indexDat<-NULL

  if(!VarCalc %in% c("None","Simulate","DeltaMethod")) VarCalc="None"

  #Subtract first year if numeric to improve convergence
  if(is.numeric(obsdat$Year) & "Year" %in% allVarNames) {
    startYear<-min(obsdat$Year)
    obsdat$Year<-obsdat$Year-startYear
    if(EstimateBycatch) logdat$Year<-logdat$Year-startYear
    if(EstimateIndex)  indexDat$Year<-indexDat$Year-startYear
  } else startYear<-min(as.numeric(as.character(obsdat$Year)))

  #Setup directory naming
  dirname<-list()
  outDir<-paste0(baseDir, paste("/Output", runName))
  if(!dir.exists(outDir)) dir.create(outDir)
  for(run in 1:numSp) {
    dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/")
    if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]])
  }

  #Make R objects to store analysis
  modelTable<-list()
  modelSelectTable<-list()
  modFits<-list()
  modPredVals<-list()
  modIndexVals<-list()
  residualTab<-list()
  bestmod<-NULL
  predbestmod<-list()
  indexbestmod<-list()
  allmods<-list()
  allindex<-list()
  modelFail<-matrix("-",numSp,length(modelTry),dimnames=list(common,modelTry))
  rmsetab<-list()
  metab<-list()





  for(run in 1:numSp) {
    residualTab[[run]]<-matrix(0,8,length(modelTry),dimnames=list(c("KS.D","KS.p", "Dispersion.ratio","Dispersion.p" ,
                                                                    "ZeroInf.ratio" ,"ZeroInf.p","Outlier" , "Outlier.p"),
                                                                  modelTry))
    modelTable[[run]]<-data.frame(model=modelTry,
                                  formula=rep("",length(modelTry)),
                                  RMSE=rep(NA,length(modelTry)),
                                  ME=rep(NA,length(modelTry)))
    modPredVals[[run]]<-rep(list(NULL),length(modelTry))
    names(modPredVals[[run]])<-modelTry
    modIndexVals[[run]]<- modPredVals[[run]]
    modFits[[run]]<- modPredVals[[run]]
    modelSelectTable[[run]]<- modPredVals[[run]]
  }





















  # list of outputs to be saved in rds file











}

