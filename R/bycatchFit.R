
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
#' @param CIval Specify confidence interval for total bycatch estimates. Should be the alpha level, e.g. 0.05 for 95%
#' @param VarCalc Character. Options are: "Simulate","DeltaMethod", or "None". Variance calculation method. Simulate will not work with a large number of sample units in the logbook data. The delta method for variance calculation is not implemented for the delta-lognormal or delta-gamma methods.
#' @param useParallel Logical. Whether to conduct the analysis using parallel processing. Only initialized when more that two cores are available.
#' @param nSims Number of simulations used to calculate confidence intervals. Ignored if \code{VarCalc} set to "None"
#' @param baseDir Character. A directory to save output. Defaults to current working directory.
#' @param plotValidation Logical. Validation. If you have true values of the total bycatch (for example in a simulation study). Make PlotValidation true and fill out the rest of the specification.
#' @param trueVals The data set that contains the true simulated total catches by year.
#' @param trueCols The column of the true simulated catches that contains true bycatch by year
#' @param reportType Character. Choose type of report with results to be produced. Options are pdf, html (default) or both.
#' @import MuMIn ggplot2 parallel dplyr doParallel foreach utils tidyverse parallelly
#' @export
#' @keywords Fitting functions
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
#' #-------------------------------------------------
#' #Step 1. Run the datasetup function and review data inputs
#'setupObj<-bycatchSetup_new(
#' obsdat = obsdatExample,
#' logdat = logdatExample,
#' yearVar = "Year",
#' obsEffort = "sampled.sets",
#' logEffort = "sets",
#' obsCatch = "Catch",
#' catchUnit = "number",
#' catchType = "dead discard"
#' logNum = NA,
#' sampleUnit = "trips",
#' factorVariables = c("Year","season"),
#' numericVariables = NA,
#' includeObsCatch  = FALSE,
#' logUnsampledEffort = NULL,
#' matchColumn = NA,
#' EstimateIndex = FALSE,
#' EstimateBycatch = TRUE,
#' baseDir = getwd(),
#' runName = "SimulatedExample",
#' runDescription = "Example with simulated data",
#' common = "Simulated species",
#' sp = "Genus species",
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
#' CIval = 0.05,
#' VarCalc = "Simulate",
#' useParallel = TRUE,
#' nSims = 1000,
#' baseDir = getwd(),
#' plotValidation = FALSE,
#' trueVals = NULL,
#' trueCols = NULL
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
  CIval = 0.05,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 10,
  baseDir = getwd(),
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL,
  reportType = "html"
){

  #unpack setup obj
  obsdat<-logdat<-yearVar<-obsEffort<-logEffort<-obsCatch<-catchUnit<-catchType<-
    logNum<-sampleUnit<-factorVariables<-numericVariables<-
    logUnsampledEffort<-includeObsCatch<-matchColumn<-
    EstimateIndex<-EstimateBycatch<-
    baseDir<-runName<-runDescription<-common<-sp<-NULL

  # some of these will be rewritten along the code
  dat<-numSp<-yearSum<-allVarNames<-startYear<-strataSum<-NULL

  for(r in 1:NROW(setupObj$bycatchInputs)) assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]])
  for(r in 1:NROW(setupObj$bycatchOutputs)) assign(names(setupObj$bycatchOutputs)[r],setupObj$bycatchOutputs[[r]])

  # #unpack designObj - if there is a designObj
  # if(!is.null(designObj)){
  # designMethods<-designVars<-designPooling<-poolTypes<-pooledVar<-adjacentNum<-minStrataUnit<-baseDir<-NULL
  #
  # yearSum<-yearSumGraph<-strataSum<-poolingSum<-includePool<-designyeardf<-designstratadf<-NULL
  #
  # for(r in 1:NROW(designObj$designInputs)) assign(names(designObj$designInputs)[r], designObj$designInputs[[r]])
  # for(r in 1:NROW(designObj$designOutputs)) assign(names(designObj$designOutputs)[r],designObj$designOutputs[[r]])
  # }

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
    dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/","bycatchFit files/")
    if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]],recursive = TRUE)
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

  ############## This is the main analysis loop ##########################
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

    datval<-dat[[run]]
    outVal<-dirname[[run]]
    varExclude<-NULL

    #Fit all models except delta
    for(mod in which(!grepl("delta",modelTry,ignore.case=TRUE))){
      modFit<-suppressWarnings(BycatchEstimator:::findBestModelFunc(
        obsdatval = datval,
        modType = modelTry[mod],
        printOutput=TRUE,
        requiredVarNames = requiredVarNames,
        allVarNames = allVarNames,
        complexModel = complexModel,
        common = common,
        randomEffects = randomEffects,
        useParallel = useParallel,
        selectCriteria = selectCriteria,
        catchType = catchType,
        varExclude = varExclude,
        dirname = dirname,
        run = run
      ))
      modelSelectTable[[run]][[modelTry[mod]]]<-modFit[[2]]
      modFits[[run]][[modelTry[mod]]]<-modFit[[1]]
    }

    #Fit delta models
    if(any(grepl("delta",modelTry,ignore.case=TRUE))) {  #Delta models if requested
      posdat<-filter(dat[[run]], .data$pres==1)
      y<-unlist(lapply(posdat[,factorVariables],function(x) length(setdiff(levels(x),x)))) #See if all levels are included
      varExclude<-names(y)[y>0]
      if(length(varExclude>0)) print(paste(common[run], "excluding variable",varExclude,"from delta models for positive catch"))
      if((min(summary(posdat$Year))>0 |  is.numeric(datval$Year)) &
         (!is.null(modFits[[run]][["Binomial"]]) | !is.null(modFits[[run]][["TMBbinomial"]]))) { #If all years have at least one positive observation and binomial converged, carry on with delta models
        for(mod in which(grepl("delta",modelTry,ignore.case=TRUE)))  {
          modFit<-suppressWarnings(BycatchEstimator:::findBestModelFunc(
            obsdatval = posdat,
            modType = modelTry[mod],
            requiredVarNames = requiredVarNames,
            allVarNames = allVarNames,
            complexModel = complexModel,
            randomEffects = randomEffects2,
            useParallel = useParallel,
            selectCriteria = selectCriteria,
            varExclude = varExclude,
            printOutput=TRUE,
            catchType = catchType,
            common = common,
            dirname = dirname,
            run = run
          ))
          modelSelectTable[[run]][[modelTry[mod]]]<-modFit[[2]]
          modFits[[run]][[modelTry[mod]]]<-modFit[[1]]
        }
      } else {
        print("Not all years have positive observations, skipping delta models")
        modelFail[run,c("Delta-Lognormal","Delta-Gamma")]<-"data"
        modPredVals[[run]][[modelTry[mod]]]<-NULL
        modIndexVals[[run]][[modelTry[mod]]]<-NULL
      }
    }

    # If estimating bycatch, see if data set is large
    if(EstimateBycatch) {
      BigData<-ifelse(sum(logdat$SampleUnits)>10000, TRUE, FALSE)
      #Add stratum designation and check sample size in strata
      if(length(requiredVarNames)>1) {
        logdat$strata<-apply( logdat[ , requiredVarNames ] , 1 , paste , collapse = "-" )
      }
      if(length(requiredVarNames)==1)   {
        logdat$strata <- pull(logdat,var=requiredVarNames)
      }
      if(length(requiredVarNames)==0)   {
        logdat$strata <- rep(1,nrow(logdat))
      }
      if(max(tapply(logdat$SampleUnits,logdat$strata,sum))>100000) {
        print("Cannot calculate variance for large number of logbook sample units")
        VarCalc<-"None"
      }
    }
    for(mod in 1:length(modelTry)) {
      if(!is.null(modFits[[run]][[modelTry[mod]]])) {
        if(grepl("delta",modelTry[mod],ignore.case = TRUE )) {
          if(grepl("TMB",modelTry[mod])) modFit1=modFits[[run]][["TMBbinomial"]] else
            modFit1=modFits[[run]][["Binomial"]]
          modFit2<-modFits[[run]][[modelTry[mod]]]
        }
        if(!grepl("delta",modelTry[mod],ignore.case = TRUE )) {
          modFit1<-modFits[[run]][[modelTry[mod]]]
          modFit2<-NULL
        }
        if(EstimateBycatch) {
          if(VarCalc=="Simulate" |(VarCalc=="DeltaMethod" & modelTry[mod] %in% c("Delta-Lognormal","Delta-Gamma","Tweedie", "TMBdelta-Lognormal","TMBdelta-Gamma")))
            modPredVals[[run]][[modelTry[mod]]]<-
              makePredictionsSimVarBig(
                modfit1=modFit1,
                modfit2=modFit2,
                modtype=modelTry[mod],
                newdat=logdat,
                obsdatval=datval,
                includeObsCatch = includeObsCatch,
                nsim = nSims,
                requiredVarNames = requiredVarNames,
                CIval = CIval,
                common = common,
                catchType = catchType,
                dirname = dirname,
                run = run,
                randomEffects=randomEffects,
                randomEffects2=randomEffects2)
          if(VarCalc=="DeltaMethod" & !modelTry[mod] %in% c("Delta-Lognormal","Delta-Gamma","Tweedie","TMBdelta-Lognormal","TMBdelta-Gamma"))
            modPredVals[[run]][[modelTry[mod]]]<-
              makePredictionsDeltaVar(
                modfit1=modFit1,
                modtype=modelTry[mod],
                newdat=logdat,
                obsdatval=datval,
                includeObsCatch = includeObsCatch,
                requiredVarNames = requiredVarNames,
                CIval = CIval,
                common = common,
                catchType = catchType,
                dirname = dirname,
                run = run
              )
          if(VarCalc=="None") {
            modPredVals[[run]][[modelTry[mod]]]<-
              makePredictionsNoVar(
                modfit1=modFit1,
                modfit2=modFit2,
                modtype=modelTry[mod],
                newdat=logdat,
                obsdatval=datval,
                includeObsCatch = includeObsCatch,
                nsims = nSims,
                requiredVarNames = requiredVarNames,
                common = common,
                catchType = catchType,
                dirname = dirname,
                run = run
              )
          }
        }
        if(EstimateIndex) {
          modIndexVals[[run]][[modelTry[mod]]]<-
            makeIndexVar(
              modfit1=modFit1,
              modfit2=modFit2,
              modType=modelTry[mod],
              indexVarNames = indexVarNames,
              newdat = indexDat,
              printOutput=TRUE,
              nsims = nSims,
              common = common,
              catchType = catchType,
              dirname = dirname,
              run = run
            )
        }
        modelTable[[run]]$formula[mod]<-paste(formula(modFits[[run]][[modelTry[mod]]]))[[3]]
        #exclude code to print pdfs - turn fileName NULL?
        #temp<-ResidualsFunc(modFits[[run]][[modelTry[mod]]],modelTry[mod],paste0(outVal,"Residuals",modelTry[mod],".pdf"))
        temp<-ResidualsFunc(modFits[[run]][[modelTry[mod]]],modelTry[mod],fileName = NULL) #no pdfs produced
        if(!is.null(temp)) {
          residualTab[[run]][,modelTry[mod]]<-temp
          #if(residualTab[[run]]["KS.p",modelTry[mod]]<0.01 & ResidualTest) modelFail[run,modelTry[mod]]<-"resid"
          # or exclude next line of code?
          # if(residualTab[[run]]["KS.p",modelTry[mod]]<0.01) modelFail[run,modelTry[mod]]<-"resid"
        }
        if(is.null(modPredVals[[run]][[modelTry[mod]]]) & EstimateBycatch) modelFail[run,modelTry[mod]]<-"cv"
      } else {
        if(modelFail[run,modelTry[mod]]=="-") modelFail[run,modelTry[mod]]<-"fit"
      }

      #if(length(which(modelFail[run,]=="-"))<1){
      #warning("No models were able to converge")
      #}

    }

    #Combine all predictions, except Binomial
    if(EstimateBycatch) {
      if(is.factor(modPredVals[[run]][[1]]$Year))
      #   yearSumGraph[[run]]$Year<-factor(yearSumGraph[[run]]$Year)
      # allmods[[run]]<-bind_rows(modPredVals[[run]],.id="Source") %>%
      #   filter(!.data$Source=="Binomial",!.data$Source=="TMBbinomial")
      # allmods[[run]]<-bind_rows(allmods[[run]],yearSumGraph[[run]])
      # allmods[[run]]$Valid<-ifelse(modelFail[run,match(allmods[[run]]$Source,dimnames(modelFail)[[2]])]=="-" | allmods[[run]]$Source %in% c("Unstratified ratio","Ratio","Design Delta"),1,0)

      allmods[[run]]<-bind_rows(modPredVals[[run]],.id="Source") %>%
        filter(!.data$Source=="Binomial",!.data$Source=="TMBbinomial")
      allmods[[run]]$Valid<-ifelse(modelFail[run,match(allmods[[run]]$Source,dimnames(modelFail)[[2]])]=="-",1,0)
    }

    if(EstimateIndex) {
      allindex[[run]]<-bind_rows(modIndexVals[[run]],.id="Source") %>%
        filter(!.data$Source=="Binomial",!.data$Source=="TMBbinonomial")
      allindex[[run]]$Valid<-ifelse(modelFail[run,match(allindex[[run]]$Source,dimnames(modelFail)[[2]])]=="-" ,1,0)
    }

    #Print the diagnostic tables
    #write.csv(residualTab[[run]],paste0(outVal,"residualDiagnostics.csv")) #being saved also as modelSummary.csv
    #write.csv(modelFail,paste0(outDir,"/modelFail.csv")) #exclude this output

    ######## Cross validation 10 fold  ####################################
    rmsetab[[run]]<-matrix(NA,10,length(modelTry),dimnames=list(1:10,modelTry))
    rmsetab[[run]]<-rmsetab[[run]][,colnames(rmsetab[[run]])!="Binomial"]
    metab[[run]]<-rmsetab[[run]]
    if(DoCrossValidation & length(which(modelFail[run,!colnames(modelFail)%in%c("Binomial","TMBbinomial")]=="-"))>1) {  #Don't do unless at least one model worked
      datval$cvsample<-sample(rep(1:10,length=dim(datval)[1]),replace=FALSE)
      for(i in 1:10 ) {
        datin<-datval[datval$cvsample!=i,]
        datout<-datval[datval$cvsample==i,]
        datout$SampleUnits<-rep(1,dim(datout)[1])
        for(mod in which(!grepl("delta",modelTry,ignore.case = TRUE))) {
          if(modelFail[run,modelTry[mod]]=="-") {
            # if(DredgeCrossValidation) { #exclude DredgeCrossValidation, but need to understand code first
            #   modFit1<-suppressWarnings(findBestModelFunc(
            #     datin,
            #     modelTry[mod],
            #     requiredVarNames = requiredVarNames,
            #     allVarNames = allVarNames,
            #     complexModel = complexModel,
            #     useParallel = useParallel,
            #     selectCriteria = selectCriteria,
            #     catchType = catchType,
            #     varExclude = varExclude,
            #     randomEffects=randomEffects
            #   ))[[1]]
            #} else {
              modFit1<-suppressWarnings(FitModelFuncCV(formula(paste0("y~",modelTable[[run]]$formula[mod])),
                                                       modType=modelTry[mod],obsdatval=datin))
            #}

            if(!modelTry[mod] %in% c("Binomial","TMBbinomial")) {
              predcpue<-makePredictions(
                modFit1,
                modType=modelTry[mod],
                newdat = datout
              )
              rmsetab[[run]][i,modelTry[mod]]<-getRMSE(predcpue$est.cpue,datout$cpue)
              metab[[run]][i,modelTry[mod]]<-getME(predcpue$est.cpue,datout$cpue)
            } else {
              if(modelTry[mod]=="Binomial")
                binomial1<-modFit1
              if(modelTry[mod]=="TMBbinomial")
                tmbbin1<-modFit1
            }
          }
        }
        if(any(grepl("delta",modelTry,ignore.case = TRUE))) {
          posdat<-filter(datin, .data$pres==1)
          y<-unlist(lapply(posdat[,factorVariables],function(x) length(setdiff(levels(x),x)))) #See if all levels are included
          varExcludecv<-names(y)[y>0]
          for(mod in which(grepl("delta",modelTry,ignore.case = TRUE))) {
            if(modelFail[run,modelTry[mod]]=="-" & !(!is.numeric(posdat$Year) & min(table(posdat$Year))==0)) {
              # if(DredgeCrossValidation) {
              #   modFit1<-suppressWarnings(findBestModelFunc(
              #     posdat,
              #     modelTry[mod],
              #     requiredVarNames = requiredVarNames,
              #     allVarNames = allVarNames,
              #     complexModel = complexModel,
              #     useParallel = useParallel,
              #     selectCriteria = selectCriteria,
              #     catchType = catchType,
              #     varExclude = varExcludecv,
              #     randomEffects=randomEffects2
              #   ))[[1]]
              #} else {
                if(length(varExcludecv)==0)
                  modFit1<-suppressWarnings(FitModelFuncCV(formula(paste0("y~",modelTable[[run]]$formula[mod])),modType=modelTry[mod],obsdatval=posdat)) else {
                    temp<-strsplit(gsub(" ","",modelTable[[run]]$formula[mod]),"[+]")[[1]]
                    temp<-paste(temp[!temp %in%  varExcludecv],collapse="+")
                    modFit1<-suppressWarnings(FitModelFuncCV(formula(paste0("y~",temp)),modType=modelTry[mod],obsdatval=posdat))
                  }
              #}

              if(grepl("TMB",modelTry[mod]))
                bin1<-tmbbin1 else
                  bin1<-binomial1
                predcpue<-makePredictions(bin1,modFit1,modelTry[mod],datout)
                rmsetab[[run]][i,modelTry[mod]]<-getRMSE(predcpue$est.cpue,datout$cpue)
                metab[[run]][i,modelTry[mod]]<-getME(predcpue$est.cpue,datout$cpue)
            }
          }
        }
      } #close loop for each fold of the cross validation

      # Calculate RMSE and ME
      modelTable[[run]]$RMSE[modelTable[[run]]$model!="Binomial"]<-apply(rmsetab[[run]],2,mean,na.rm=TRUE)
      modelTable[[run]]$ME[modelTable[[run]]$model!="Binomial"]<-apply(metab[[run]],2,mean,na.rm=TRUE)
      write.csv(residualTab[[run]],paste0(outVal,"modelSummary.csv"))
      write.csv(rmsetab[[run]],paste0(outVal,"rmse.csv"))
      write.csv(metab[[run]],paste0(outVal,"me.csv"))
      write.csv(modelTable[[run]],paste0(outVal,"crossvalSummary.csv"))
      #Select best model based on cross validation
      best<-which(!is.na( modelTable[[run]]$RMSE) &
                    modelTable[[run]]$RMSE==min(modelTable[[run]]$RMSE,na.rm=TRUE))
      if(length(best) >0 ) {
        bestmod[run]<-modelTry[best]
        predbestmod[[run]]<-modPredVals[[run]][[modelTry[best]]]
        indexbestmod[[run]]<-modIndexVals[[run]][[modelTry[best]]]
      } else {
        bestmod[run]<-"None"
      }
    }

    print(paste(run, common[run],"complete, ",Sys.time()))

  } #close loop for each spp

  # list of inputs/outputs to be saved in rds file
  output <- list(
    modelInputs = list(
      modelTry = modelTry,
      randomEffects = randomEffects,
      randomEffects2 = randomEffects2,
      complexModel = complexModel,
      simpleModel = simpleModel,
      indexModel = indexModel,
      selectCriteria = selectCriteria,
      DoCrossValidation = DoCrossValidation,
      CIval = CIval,
      VarCalc = VarCalc,
      useParallel = useParallel,
      nSims = nSims,
      plotValidation = plotValidation,
      trueVals = trueVals,
      trueCols = trueCols
    ),
    modelOutputs = list(
      modelTable = modelTable,
      modelSelectTable = modelSelectTable,
      modFits = modFits,
      modPredVals = modPredVals,
      modIndexVals = modIndexVals,
      indexDat = indexDat,
      indexVarNames = indexVarNames,
      residualTab = residualTab,
      bestmod = bestmod,
      predbestmod = predbestmod,
      indexbestmod = indexbestmod,
      allmods = allmods,
      allindex = allindex,
      modelFail = modelFail,
      rmsetab = rmsetab,
      metab = metab,
      dat = dat,
      requiredVarNames = requiredVarNames,
      allVarNames = allVarNames,
      startYear = startYear,
      NumCores = NumCores
    )
  )

  saveRDS(output, file=paste0(outDir,"/", Sys.Date(),"_BycatchFitSpecification.rds"))

  #Create report
  for(run in 1:numSp) {
    if(reportType == "html" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "PrintBycatchFit.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )

      if(!is.null( mkd)){

        rmarkdown::render(mkd,
                          params=list(outDir=outDir, run = run),
                          output_format = "html_document",
                          output_file = paste0(common[run], " ",catchType[run], " Model-based estimation results.html"),
                          output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                          quiet = TRUE)
      }
    }

    if(reportType == "pdf" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "PrintBycatchFit.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )

      if(!is.null( mkd)){

        tryCatch({
          rmarkdown::render(mkd,
                            params=list(outDir=outDir, run = run),
                            output_format = "pdf_document",
                            output_file = paste0(common[run], " ",catchType[run], " Model-based estimation results.pdf"),
                            output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                            quiet = TRUE)

        },
        error = function(e){
          message("PDF rendering failed, reverting to html.")
          rmarkdown::render(mkd,
                            params=list(outDir=outDir, run = run),
                            output_format = "html_document",
                            output_file = paste0(common[run], " ",catchType[run], " Model-based estimation results.html"),
                            output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                            quiet = TRUE)
        })
      }
    }
  }

  #Clean up: delete the figures/ directory after rendering
  fig_dir <- file.path(dirname(mkd), "figures")
  if (dir.exists(fig_dir)) {
    unlink(fig_dir, recursive = TRUE)
  }

  return(output)


} #close loop bycatchFit function

