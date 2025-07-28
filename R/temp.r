
PredictionsByGroup<-function(outputObj,
                             modelScenarios=NULL,
                             modelTypes=NULL,
                             spNums=NULL,
                             VarCalc=NULL,
                             groupVars)  {
  tempRunName<-abbreviate(outputObj$runName,minlength=10)
  outDir<-paste0(outputObj$baseDir, paste0("/Output", tempRunName))
  if(!dir.exists(outDir)) stop(paste("Directory",outDir,"not found."))
  modelTry<- c("Tweedie","Lognormal","Delta-Lognormal","Delta-Gamma", "TMBnbinom1","TMBlognormal",
               "TMBnbinom2","TMBtweedie","Normal","Binomial","NegBin", "TMBgamma","Gamma",
               "TMBbinomial","TMBnormal","TMBdelta-Lognormal","TMBdelta-Gamma","Poisson","TMBpoisson")
  if(is.null(modelScenarios))
    modelScenarios<-unique(outputObj$allYearEstimates$Scenario[outputObj$allYearEstimates$Source %in% modelTry])
  temp<-modelScenarios[!modelScenarios %in% outputObj$allYearEstimates$Scenario]
  if(length(temp)>0) stop(paste("ModelScenario",temp,"not found."))
  if(is.null(modelTypes))
    modelTypes=modelTry[modelTry %in% outputObj$allYearEstimates$Source]
  temp<-modelTypes[!modelTypes %in% modelTry]
  if(length(temp)>0) stop(paste("Models",temp,"not recognized."))
  temp<-modelTypes[!modelTypes %in% outputObj$allYearEstimates$Source]
  if(length(temp)>0) stop(paste("Models",temp,"not found in outputs."))
  if(is.null(spNums)) spNums=1:length(unique(outputObj$allYearEstimates$spNum))
  temp<-spNums[!spNums %in% outputObj$allYearEstimates$spNum]
  if(length(temp)>0) stop(paste("spNum",temp,"not found in outputs."))
  modTab<-expand.grid(spNum=spNums,Scenario=modelScenarios,modelType=modelTypes)
  for(i in 1:nrow(temp)) {
     setup1<-outputObj$setupObj$bycatchInputs
     run<-modTab$spNum[i]
     dirname<-gsub("Setup files","Fit files",setupObj$bycatchInputs$dirname[run])
     dirname<-paste0(dirname,abbreviate(paste(c(VarCalc,groupVars),collapse=""),minlength=20))
     if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]],recursive = TRUE)
     modelTry<-modTab$modelType[i]
     modObj1<-outputObj$modelobjList[[modTab$Scenario[i]]]
     logdat<-modObj1$modelOutputs$logdatFit
     temp<-groupVars[!groupVars %in% names(logdat)]
     if(length(temp)>0) stop("Variable",temp,"not found in logdat")
     modFits<-modObj1$modelOutputs$modFits[[run]]
     datval<-modObj1$modelOutputs$dat[[run]]
     if(is.null(VarCalc)) VarCalc<-modObj1$modelInputs$VarCalc
     if(!VarCalc %in% c("DeltaMethod","Simulate","None")) stop("VarCalc method not found.")
     if(!is.null(modFits[[modelTry]])) {
         if(grepl("delta",modelTry,ignore.case = TRUE )) {
           if(grepl("TMB",modelTry)) modFit1=modFits[["TMBbinomial"]] else
             modFit1=modFits[["Binomial"]]
           modFit2<-modFits[[modelTry]]
         }
      if(!grepl("delta",modelTry,ignore.case = TRUE )) {
           modFit1<-modFits[[modelTry]]
           modFit2<-NULL
         }
      if(VarCalc=="Simulate" |(VarCalc=="DeltaMethod" & modelTry %in% c("Delta-Lognormal","Delta-Gamma","Tweedie", "TMBdelta-Lognormal","TMBdelta-Gamma")))
             predVals<-
               makePredictionsSimVarBig(
                 modfit1=modFit1,
                 modfit2=modFit2,
                 modtype=modelTry,
                 newdat=logdat,
                 obsdatval=datval,
                 includeObsCatch =  modObj1$modelInputs$includeObsCatch,
                 nsim = modObj1$modelInputs$nSims,
                 requiredVarNames = groupVars,
                 CIval = modObj1$modelInputs$CIval,
                 common = modObj1$modelInputs$common[run],
                 catchType = modObj1$modelInputs$catchType[run],
                 dirname = dirname,
                 shortName = modObj1$modelInputs$shortName[run],
                 run = run,
                 randomEffects=modObj1$modelInputs$randomEffects,
                 randomEffects2=modObj1$modelInputs$randomEffects2,
                 modelScenario=modObj1$modelInputs$modelScenario,
                 startYear=modObj1$modelInputs$startYear)
       if(VarCalc=="DeltaMethod" & !modelTry %in% c("Delta-Lognormal","Delta-Gamma","Tweedie","TMBdelta-Lognormal","TMBdelta-Gamma"))
             predVals<-
               makePredictionsDeltaVar(
                 modfit1=modFit1,
                 modtype=modelTry,
                 newdat=logdat,
                 obsdatval=datval,
                 includeObsCatch = includeObsCatch,
                 requiredVarNames = requiredVarNames,
                 CIval = CIval,
                 common = common,
                 shortName = shortName,
                 catchType = catchType,
                 dirname = dirname,
                 run = run,
                 modelScenario=modelScenario,
                 startYear=startYear
               )
           if(VarCalc=="None") {
             modPredVals[[run]][[modelTry]]<-
               makePredictionsNoVar(
                 modfit1=modFit1,
                 modfit2=modFit2,
                 modtype=modelTry,
                 newdat=logdat,
                 obsdatval=datval,
                 includeObsCatch = includeObsCatch,
                 nsims = nSims,
                 requiredVarNames = requiredVarNames,
                 common = common,
                 catchType = catchType,
                 shortName = shortName,
                 dirname = dirname,
                 run = run,
                 modelScenario=modelScenario,
                 startYear=startYear
               )
           }






  }
}


PredictionsByGroup(outputObj,
                             modelScenarios=NULL,
                             modelTypes=NULL,
                             spNums=1)
