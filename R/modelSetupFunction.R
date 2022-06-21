

#---------------------------------
#Model setup
#----------------------------------

#Roxygen header
#'Bycatch estimation model setup
#'
#'Sets global conditions and makes a preliminary data summary.
#'
#'
#' @param modelTry  Specify which observation error models to try. Options are: "Binomial", "Normal","Lognormal", "Delta-Lognormal", "Delta-Gamma", "NegBin" for Negative binomial" using glm.mb in the MASS library, "Tweedie" for Tweedie GLM from the cpglm library, and "TMBnbinom1", "TMBnbinom2", and "TMBtweedie" for negative binomial 1, negative binomial 2 and Tweedie from the GLMMTMB library. Binomial is run automatically as part of the delta models if either of them are selected.
#' @param obsdat Observer data set
#' @param logdat Logbook data set
#' @param yearVar Character. The name of the year variable in \code{obsdat} and \code{logdat}. Both input files must contain the same variable name for year.
#' @param obsEffort Character. The name of the effort variable in \code{obsdat}. This variable must have the same effort units as \code{logEffort}
#' @param logEffort Character. The name of the effort variable in \code{logdat}. Optional and only used when estimating bycatch. This variable must have the same effort units as \code{obsEffort}
#' @param logUnsampledEffort Character. The name of the unsampled effort variable in in \code{logdat}. Optional and used to specify a column for effort that is not sampled, in trips with observers. This can be zero in all cases if observers sample 100% of effort in sampled trips. Only used when \code{includeObsCatch} is TRUE
#' @param includeObsCatch Logical. Set to TRUE if (1) the observed sample units can be matched to the logbook sample units and (2) you want to calculate total bycatch as the observed bycatch plus the predicted unobserved bycatch. This doesn't work with aggregated logbook effort.
#' @param matchColumn Character. If \code{includeObsCatch} is TRUE, give the name of the column that matches sample units between the observer and logbook data. Otherwise, this can be NA
#' @param factorNames Character vector. Specify which variables should be interpreted as categorical, ensuring imposes factor format on these variables. Variables not in this list will retain their original format. These variables must have identical names and factor levels in \code{obsdat} and \code{logdat}
#' @param logNum Character vector. The name of the column in \code{logdat} that gives the number of sample units (e.g., trips or sets). If the logbook data is not aggregated (i.e. each row is a sample unit) set value to NA
#' @param sampleUnit Character. What is the sample unit in \code{logdat}? e.g. sets or trips.
#' @param EstimateIndex Logical. What would you like to estimate? You may calculate either an annual abundance index, or total bycatch, or both.
#' @param EstimateBycatch Logical. What would you like to estimate? You may calculate either an annual abundance index, or total bycatch, or both. If you want total bycatch, you must provide logbook data or some other source of total effort to \code{logdat}.
#' @param complexModel Specify as stats::formula. Specify the most complex and simplest model to be considered. The code will find compare all intermediate models using information criteria.
#' @param simpleModel Specify as stats::formula.
#' @param indexModel Specify as stats::formula. Use indexModel to specify which strata to keep separate in calculating abundance indices.
#' @param baseDir Character. A directory to save output. Defaults to current working directory.
#' @param runName Characer. Give a name to the run, which will be used to set up a directory for the outputs
#' @param runDescription Character. Brief summary of the run, which will be used to set up a directory for the outputs
#' @param common Character vector. Provide a common name for the species used in bycatch and index estimation. Can be a vector of names to do multiple species at the same time.
#' @param sp  Character vector. Provide a scientific name for the species used in bycatch and index estimation. Can be a vector of names to do multiple species at the same time
#' @param obsCatch Character vector. The name of the column(s) in \code{obsdat} that contain catch. If it is a vector, order of variable names must follow the same order as names provided in \code{common} and \code{sp}
#' @param catchUnit Character vector. Give units of catch (e.g., number) to go in plot labels. Must be a vector of the same length as \code{sp}
#' @param catchType Character vector. Give type of catch (e.g., dead discards) to go in plot labels. Must be a vector of the same length as \code{sp}
#' @import ggplot2 parallel dplyr doParallel foreach utils kableExtra tidyverse
#' @importFrom stats median
#' @export
#' @keywords Fitting functions
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
#' setupObj<-bycatchSetup(
#' modelTry = c("Delta-Lognormal","Delta-Gamma","TMBnbinom1","TMBnbinom2","TMBtweedie"),
#' obsdat = obsdatExample,
#' logdat = logdatExample,
#' yearVar = "Year",
#' obsEffort = "sampled.sets",
#' logEffort = "sets",
#' logUnsampledEffort = NULL,
#' includeObsCatch  = FALSE,
#' matchColumn = NA,
#' factorNames = c("Year","season"),
#' EstimateIndex = TRUE,
#' EstimateBycatch = TRUE,
#' logNum = NA,
#' sampleUnit = "trips",
#' complexModel = formula(y~(Year+season)^2),
#' simpleModel = formula(y~Year),
#' indexModel = formula(y~Year),
#' baseDir = getwd(),
#' runName = "SimulatedExample",
#' runDescription = "Example with simulated data",
#' common = "Simulated species",
#' sp = "Genus species",
#' obsCatch = "Catch",
#' catchUnit = "number",
#' catchType = "dead discard"
#')}

bycatchSetup <- function(
  modelTry = c("Delta-Lognormal","Delta-Gamma","TMBnbinom1","TMBnbinom2","TMBtweedie"),
  obsdat,
  logdat,
  yearVar,
  obsEffort,
  logEffort,
  logUnsampledEffort = NULL,
  includeObsCatch  = FALSE,
  matchColumn = NA,
  factorNames,
  EstimateIndex,
  EstimateBycatch,
  logNum,
  sampleUnit,
  complexModel,
  simpleModel,
  indexModel,
  baseDir = getwd(),
  runName,
  runDescription,
  common,
  sp,
  obsCatch,
  catchUnit,
  catchType


  ){

  SampleUnits<-Year<-drop_na<-Catch<-Effort<-cpue<-pres<-Pos<-OUnit<-OEff<-Eff<-Units<-outDir<-NULL

  #Set global conditions
  theme_set(theme_bw()) #ggplot theme
  defaultW <- 0
  options(warn=defaultW)
  NumCores<-detectCores()  #Check if machine has multiple cores for parallel processing

  #Check that all models in modelTry are valid
  if(!all(modelTry %in% c("Tweedie","Lognormal","Delta-Lognormal","Delta-Gamma","TMBnbinom1","TMBnbinom2","TMBtweedie","Normal","Binomial","NegBin") ))
    stop(paste("Model requested in modelTry not available"))

   #Make sure binomial is included if either of the delta models is
  if(("Delta-Lognormal" %in% modelTry |"Delta-Gamma" %in% modelTry) & !"Binomial" %in% modelTry)
    modelTry<-c("Binomial",modelTry)

  # Set up variables
  #if("Year" %in% names(obsdat) & !yearVar=="Year") obsdat<-obsdat %>% rename(oldYear=Year)
  #if("Year" %in% names(logdat) & !yearVar=="Year") logdat<-logdat %>% rename(oldYear=Year)
  obsdat<-obsdat %>% ungroup() %>%
    rename(Year=!!yearVar)
  logdat<-logdat %>%  ungroup() %>%
    rename(Year=!!yearVar)
  requiredVarNames<-as.vector(getAllTerms(simpleModel))
  allVarNames<-as.vector(getAllTerms(complexModel))
  allVarNames<-allVarNames[grep(":",allVarNames,invert=TRUE)]
  allVarNames<-allVarNames[grep("I(*)",allVarNames,invert=TRUE)]
  if(!all(allVarNames %in% names(obsdat)))
    print(paste0("Variable ", allVarNames[!allVarNames%in% names(obsdat) ], " not found in observer data"))
  if(!all(allVarNames %in% names(logdat)))
    print(paste0("Variable ", allVarNames[!allVarNames%in% names(logdat) ], " not found in logbook data"))
  #It's all right not to see the variable name if it is a function of another variable that is present
  indexVarNames<-as.vector(getAllTerms(indexModel))
  if(!"Year" %in% indexVarNames) indexVarNames<-c("Year",indexVarNames)

  #Set up data frames
  obsdat<-obsdat %>%
    rename(Effort=!!obsEffort) %>%
    mutate_at(vars(all_of(factorNames)),factor)
  if(EstimateBycatch) {
    if(is.na(logNum))   {
      logdat<-mutate(logdat,SampleUnits=1)
      logNum<-"SampleUnits"
    }
    logdat<-logdat %>%
      rename(Effort=!!logEffort,SampleUnits=!!logNum) %>%
      mutate_at(vars(all_of(factorNames)),factor)
    if(logEffort==sampleUnit) logdat<-mutate(logdat,Effort=SampleUnits)
    if(includeObsCatch) {
      obsdat<-obsdat %>% rename(matchColumn=!!matchColumn)
      logdat<-logdat %>% rename(matchColumn=!!matchColumn,unsampledEffort=!!logUnsampledEffort)
    }
  }

  #Add stratum designation and check sample size in strata
  if(length(requiredVarNames) > 1) {
    logdat$strata<-apply( logdat[ , requiredVarNames ] , 1 , paste , collapse = "-" )
  } else {
    logdat$strata <- pull(logdat,var=requiredVarNames)
  }
  if(max(tapply(logdat$SampleUnits,logdat$strata,sum)) > 100000) {
    print("Unadvisable to calculate variance of estimators for such large number of logbook sample units")
  }

  #indexDat for making index
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

  #Subtract first year if numeric to improve convergence
  if(is.numeric(obsdat$Year) & "Year" %in% allVarNames) {
    startYear<-min(obsdat$Year)
    obsdat$Year<-obsdat$Year-startYear
    logdat$Year<-logdat$Year-startYear
    indexDat$Year<-indexDat$Year-startYear
  } else startYear<-min(as.numeric(as.character(obsdat$Year)))

  #Set up directory for output
  numSp<-length(sp)
  outDir<-paste0(baseDir, paste("/Output", runName))
  if(!dir.exists(outDir)) dir.create(outDir)

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
  strataSum<-list()

  #Make lists to keep output, which will also be output as .pdf and .csv files for use in reports.
  dirname<-list()
  dat<-list()
  #Loop through all species and print data summary. Note that records with NA in either catch or effort are excluded automatically
  yearSum<-list()
  for(run in 1:numSp) {
    dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/")
    if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]])
    if(includeObsCatch) tempvars<-c(allVarNames,"Effort","Catch","matchColumn") else
      tempvars<-c(allVarNames,"Effort","Catch")
    if(!"Year" %in%  tempvars) tempvars<-c("Year",tempvars)
    dat[[run]]<-obsdat %>%
      rename(Catch=!!obsCatch[run])%>%
      dplyr::select_at(all_of(tempvars)) %>%
      drop_na()   %>%
      mutate(cpue=Catch/Effort,
             log.cpue=log(Catch/Effort),
             pres=ifelse(cpue>0,1,0))
    if(dim(dat[[run]])[1]<dim(obsdat)[1]) print(paste0("Removed ",dim(obsdat)[1]-dim(dat[[run]])[1]," rows with NA values for ",common[run]))
    # yearSum[[run]]<-dat[[run]] %>% group_by(Year) %>%
    #   summarize(OCat=sum(Catch,na.rm=TRUE),
    #             OEff=sum(Effort,na.rm=TRUE),
    #             OUnit=length(Year),
    #             CPUE=mean(cpue,na.rm=TRUE),
    #             CPse=standard.error(cpue),
    #             Out=outlierCountFunc(cpue),
    #             Pos=sum(pres,na.rm=TRUE)) %>%
    #   mutate(PFrac=Pos/OUnit)
    # if(EstimateBycatch) {
    #   x<-logdat  %>% group_by(Year) %>%
    #     summarize(Eff=sum(Effort,na.rm=TRUE),Units=sum(SampleUnits))
    #   yearSum[[run]]<-merge(yearSum[[run]],x) %>% mutate(EFrac=OEff/Eff,
    #                                                      UFrac=OUnit/Units)
    #   logyear<-logdat %>% group_by(Year) %>% summarize(Effort=sum(Effort,na.rm=TRUE))
    #   x=ratio.func(dat[[run]]$Effort,dat[[run]]$Catch,dat[[run]]$Year,
    #                logyear$Effort,logyear$Effort,logyear$Year)
    #   yearSum[[run]]<-cbind(yearSum[[run]],Cat=x$stratum.est,Cse=x$stratum.se) %>%
    #     ungroup() %>% mutate(Year=as.numeric(as.character(Year))) %>%
    #     mutate(Year=ifelse(Year<startYear,Year+startYear,Year))
    # }
    #write.csv(yearSum[[run]],paste(dirname[[run]],common[run],catchType[run],"DataSummary.csv"))

    yearSum[[run]]<-MakeSummary(
                      obsdatval = dat[[run]],
                      logdatval = logdat,
                      strataVars = "Year",
                      EstimateBycatch = EstimateBycatch,
                      startYear = startYear
                    )
    write.csv(dplyr::select(yearSum[[run]],c("Year","OCat","OEff","OUnit","CPUE","CPse","Out","Pos","PFrac","Eff","Units","EFrac","UFrac","Cat","Cse")),
              paste0(dirname[[run]],common[run],catchType[run],"DataSummary.csv"))
    strataSum[[run]]<-MakeSummary(
                        obsdatval = dat[[run]],
                        logdatval = logdat,
                        strataVars = requiredVarNames,
                        EstimateBycatch = EstimateBycatch,
                        startYear = startYear
                      )
    write.csv(strataSum[[run]],
              paste0(dirname[[run]],common[run],catchType[run],"StrataSummary.csv"))
  }

  #Create report
  save(list=c("numSp","yearSum","runName", "common", "sp"),file=paste0(outDir,"/","sumdatR"))
  mkd<-tryCatch({
    system.file("Markdown", "PrintDataSummary.Rmd", package = "BycatchEstimator", mustWork = TRUE)
  },
    error = function(c) NULL
  )
  if(!is.null( mkd)){
    rmarkdown::render(mkd,
                      params=list(OutDir=outDir),
                      output_file = "DataSummary.pdf",
                      output_dir=outDir,
                      quiet = TRUE)
  }

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

  #Create output list
  output<-list(

    #Inputs to the model
    bycatchInputs = list(
      modelTry = modelTry,
      obsdat = obsdat,
      logdat = logdat,
      yearVar = yearVar,
      obsEffort = obsEffort,
      logEffort = logEffort,
      logUnsampledEffort = logUnsampledEffort,
      includeObsCatch  = includeObsCatch,
      matchColumn = matchColumn,
      factorNames = factorNames,
      EstimateIndex = EstimateIndex,
      EstimateBycatch =EstimateBycatch,
      logNum = logNum,
      sampleUnit = sampleUnit,
      complexModel = complexModel,
      simpleModel = simpleModel,
      indexModel = indexModel,
      baseDir = baseDir,
      runName = runName,
      runDescription = runDescription,
      common = common,
      sp = sp,
      obsCatch = obsCatch,
      catchUnit = catchUnit,
      catchType = catchType
    ),

    #Outputs from data processing
    bycatchOutputs = list(
      numSp = numSp,
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
      metab = rmsetab,
      dat = dat,
      numSp = numSp,
      yearSum = yearSum,
      requiredVarNames = requiredVarNames,
      allVarNames = allVarNames,
      startYear = startYear,
      strataSum = strataSum
    )
  )

  saveRDS(output, file=paste0(outDir,"/", Sys.Date(),"_BycatchModelSpecification.rds"))
  return(output)
}

















