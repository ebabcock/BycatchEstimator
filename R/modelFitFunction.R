

#---------------------------------
#Model setup
#----------------------------------

#Roxygen header
#'Bycatch model fitting routine
#'
#'Produces model-based estimates of bycatch and annual abundance index, as specified in \code{bycatchSetup}
#'
#'
#' @param setupObj  An object produced by \code{bycatchSetup}.
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
#' @import MuMIn
#' @export
#' @keywords Fitting functions
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
#' #-------------------------------------------------
#' #Step 1. Run the setup file and review data inputs
#'setupObj<-bycatchSetup(
#'  modelTry = c("Lognormal","Delta-Lognormal","Delta-Gamma","TMBnbinom1","TMBnbinom2","TMBtweedie"),
#'  obsdat = LLSIM_BUM_Example_observer,
#'  logdat = LLSIM_BUM_Example_logbook,
#'  yearVar = "Year",
#'  obsEffort = "hooks",
#'  logEffort = "hooks",
#'  logUnsampledEffort = "unsampledEffort",
#'  includeObsCatch  = TRUE,
#'  matchColumn = "trip",
#'  factorNames = c("Year","fleet","area","season"),
#'  EstimateIndex = TRUE,
#'  EstimateBycatch = TRUE,
#'  logNum = NA,
#'  sampleUnit = "trips",
#'  complexModel = formula(y~Year+fleet+hbf+area+season+Year:area),
#'  simpleModel = formula(y~Year+fleet+area),
#'  indexModel = formula(y~Year+area),
#'  baseDir = getwd(),
#'  runName = "LLSIMBUMtrip2022Aprilobs05mc",
#'  runDescription = "LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022",
#'  common = c("Swordfish","Blue marlin")[2],
#'  sp = c("Xiphias gladius","Makaira nigricans")[2],
#'  obsCatch = c("SWO","BUM")[2],
#'  catchUnit = "number",
#'  catchType = "catch"
#')
#'
#'-------------
#' #Step 2. Model Fitting
#'bycatchFit(
#'  setupObj = setupObj,
#'  selectCriteria = "BIC",
#'  DoCrossValidation = TRUE,
#'  DredgeCrossValidation = FALSE,
#'  ResidualTest = FALSE,
#'  CIval = 0.05,
#'  VarCalc = "Simulate",
#'  useParallel = TRUE,
#'  nSims = 1000,
#'  baseDir = getwd(),
#'  plotValidation = FALSE,
#'  trueVals = NULL,
#'  trueCols = NULL
#')}


bycatchFit<-function(

  setupObj,
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

  #Unpack setupObj
  modelTry<-obsdat<-logdat<-yearVar<-obsEffort<-logEffort<-logUnsampledEffort<-
    includeObsCatch<-matchColumn<-factorNames<-randomEffects<-
  EstimateIndex<-EstimateBycatch<-logNum<-sampleUnit<-complexModel<-simpleModel<-indexModel<-
    designMethods<-designVars<-designPooling<-minStrataUnit<-minStrataEffort<-
    baseDir<-runName<-runDescription<-
  common<-sp<-obsCatch<-catchUnit<-catchType<-NULL

  numSp<-modelTable<-modelSelectTable<-modFits<-modPredVals<-modIndexVals<-residualTab<-bestmod<-predbestmod<-indexbestmod<-allmods<-allindex<-
  modelFail<-rmsetab<-metab<-dat<-yearSum<-requiredVarNames<-allVarNames<-indexDat<-strataSum<-NumCores<-NULL

  for(r in 1:NROW(setupObj$bycatchInputs)) assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]])
  for(r in 1:NROW(setupObj$bycatchOutputs)) assign(names(setupObj$bycatchOutputs)[r],setupObj$bycatchOutputs[[r]])

  #Setup directory naming
  dirname<-list()
  outDir<-paste0(baseDir, paste("/Output", runName))
  if(!dir.exists(outDir)) dir.create(outDir)
  for(run in 1:numSp) {
    dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/")
    if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]])
  }

  #Make sure there are multiple cores to use Parallel processing
  if(NumCores<=1) useParallel=FALSE

  ############## This is the main analysis loop ##########################
  for(run in 1:numSp) {
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
      y<-unlist(lapply(posdat[,factorNames],function(x) length(setdiff(levels(x),x)))) #See if all levels are included
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
            randomEffects = randomEffects,
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

    # See if data set is large
    BigData<-ifelse(sum(logdat$SampleUnits)>10000, TRUE, FALSE)
    #Add stratum designation and check sample size in strata
    if(length(requiredVarNames)>1) {
      logdat$strata<-apply( logdat[ , requiredVarNames ] , 1 , paste , collapse = "-" )
    } else {
      logdat$strata <- pull(logdat,var=requiredVarNames)
    }
    if(max(tapply(logdat$SampleUnits,logdat$strata,sum))>100000) {
      print("Cannot calculate variance for large number of logbook sample units")
      VarCalc<-"None"
    }

    #Make predictions, residuals, etc. for all models
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
          if(VarCalc=="Simulate" |(VarCalc=="DeltaMethod" & modelTry[mod] %in% c("Delta-Lognormal","Delta-Gamma","Tweedie", "TMBdelta-lognormal","TMBdelta-gamma")))
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
                  run = run)
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
        temp<-ResidualsFunc(modFits[[run]][[modelTry[mod]]],modelTry[mod],paste0(outVal,"Residuals",modelTry[mod],".pdf"))
        if(!is.null(temp)) {
          residualTab[[run]][,modelTry[mod]]<-temp
          if(residualTab[[run]]["KS.p",modelTry[mod]]<0.01 & ResidualTest) modelFail[run,modelTry[mod]]<-"resid"
        }
        if(is.null(modPredVals[[run]][[modelTry[mod]]])) modelFail[run,modelTry[mod]]<-"cv"
      } else {
        if(modelFail[run,modelTry[mod]]=="-") modelFail[run,modelTry[mod]]<-"fit"
      }
    }

    #Combine all predictions, except Binomial
    if(EstimateBycatch) {
      if(is.factor(modPredVals[[run]][[1]]$Year))
        yearSumGraph[[run]]$Year<-factor(yearSumGraph[[run]]$Year)
       allmods[[run]]<-bind_rows(modPredVals[[run]],.id="Source") %>%
        filter(!.data$Source=="Binomial",!.data$Source=="TMBbinomial")
      allmods[[run]]<-bind_rows(allmods[[run]],yearSumGraph[[run]])
      allmods[[run]]$Valid<-ifelse(modelFail[run,match(allmods[[run]]$Source,dimnames(modelFail)[[2]])]=="-" | allmods[[run]]$Source %in% c("Unstratified ratio","Ratio","Design Delta"),1,0)
    }
    if(EstimateIndex) {
      allindex[[run]]<-bind_rows(modIndexVals[[run]],.id="Source") %>%
        filter(!.data$Source=="Binomial",!.data$Source=="TMBbinomial")
      allindex[[run]]$Valid<-ifelse(modelFail[run,match(allindex[[run]]$Source,dimnames(modelFail)[[2]])]=="-" ,1,0)
    }

    #Print the diagnostic tables
    write.csv(residualTab[[run]],paste0(outVal,"residualDiagnostics.csv"))
    write.csv(modelFail,paste0(outDir,"/modelFail.csv"))

    ######## Cross validation 10 fold  ####################################
    rmsetab[[run]]<-matrix(NA,10,length(modelTry),dimnames=list(1:10,modelTry))
    rmsetab[[run]]<-rmsetab[[run]][,colnames(rmsetab[[run]])!="Binomial"]
    metab[[run]]<-rmsetab[[run]]
    if(DoCrossValidation & length(which(modelFail[run,!colnames(modelFail)%in%c("Binomial","TMBbinomial")]=="-"))>1) {  #Don't do unless at least one model worked
      datval$cvsample<-sample(rep(1:10,length=dim(datval)[1]),replace=FALSE)
      table(datval$cvsample,datval$Year)
      for(i in 1:10 ) {
        datin<-datval[datval$cvsample!=i,]
        datout<-datval[datval$cvsample==i,]
        datout$SampleUnits<-rep(1,dim(datout)[1])
        for(mod in which(!grepl("delta",modelTry,ignore.case = TRUE))) {
          if(modelFail[run,modelTry[mod]]=="-") {
            if(DredgeCrossValidation) {
              modFit1<-suppressWarnings(findBestModelFunc(
                datin,
                modelTry[mod],
                requiredVarNames = requiredVarNames,
                allVarNames = allVarNames,
                complexModel = complexModel,
                randomEffects = randomEffects,
                useParallel = useParallel,
                selectCriteria = selectCriteria,
                catchType = catchType,
                varExclude = varExclude
              ))[[1]]
            } else {
              modFit1<-suppressWarnings(FitModelFuncCV(formula(paste0("y~",modelTable[[run]]$formula[mod])),
                                      modType=modelTry[mod],obsdatval=datin))
            }
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
                 tmbbin1<-modfit1
            }
          }
        }
        if(any(grepl("delta",modelTry,ignore.case = TRUE))) {
          posdat<-filter(datin, .data$pres==1)
          for(mod in which(grepl("delta",modelTry,ignore.case = TRUE))) {
            if(modelFail[run,modelTry[mod]]=="-" & !(!is.numeric(posdat$Year) & min(table(posdat$Year))==0)) {
              if(DredgeCrossValidation) {
                modFit1<-suppressWarnings(findBestModelFunc(
                  posdat,
                  modelTry[mod],
                  requiredVarNames = requiredVarNames,
                  allVarNames = allVarNames,
                  complexModel = complexModel,
                  useParallel = useParallel,
                  selectCriteria = selectCriteria,
                  catchType = catchType,
                  varExclude = varExclude
                ))[[1]]
              } else {
                modFit1<-suppressWarnings(FitModelFuncCV(formula(paste0("y~",modelTable[[run]]$formula[mod])),modType=modelTry[mod],obsdatval=posdat))
              }
              bin1<-ifelse(modelTry[mod]=="Binomial",binomial1,tmbbin1)
              predcpue<-makePredictions(bin1,modFit1,modelTry[mod],datout)
              rmsetab[[run]][i,modelTry[mod]]<-getRMSE(predcpue$est.cpue,datout$cpue)
              metab[[run]][i,modelTry[mod]]<-getME(predcpue$est.cpue,datout$cpue)
            }
          }
        }
      }
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

    if(doReport){
      save(list=c("numSp","yearSum","runName","runDescription","allVarNames",
                  "common", "sp","bestmod","CIval","includeObsCatch",
                  "predbestmod","indexbestmod","allmods","allindex","modelTable",
                  "modelSelectTable","modFits","modPredVals","VarCalc"
                  ,"modIndexVals","modelFail","rmsetab","metab",
                  "residualTab" ,"run","modelTry","EstimateIndex","EstimateBycatch",
                  "DoCrossValidation","indexVarNames","selectCriteria","sampleUnit",
                  "modelTry","catchType","catchUnit","residualTab",
                  "plotValidation","trueVals","trueCols","startYear"),
           file=paste0(outVal,"/","resultsR"))

      mkd<-tryCatch({
        system.file("Markdown", "PrintResults.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )
      if(!is.null( mkd)){
        suppressWarnings(rmarkdown::render(mkd,
                          params=list(outVal=outVal),
                          output_file = paste0(common[run], catchType[run],"results.pdf"),
                          output_dir=outVal,
                          quiet = TRUE))
      }
    }
    print(paste(run, common[run],"complete, ",Sys.time()))

  }
}
