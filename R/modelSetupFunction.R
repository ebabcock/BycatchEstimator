
#---------------------------------
#Model setup
#----------------------------------

#Roxygen header
#'Bycatch estimation model setup
#'
#'Sets global conditions and makes a preliminary data summary.
#'
#'
#' @param modelTry  Specify which observation error models to try. Options are: "Binomial", "Normal","Lognormal", "Delta-Lognormal", and "Delta-Gamma", for models using the lm and glm functions, "NegBin" for Negative binomial using glm.nb in the MASS library, "Tweedie" for Tweedie GLM from the cpglm function in the cplm library, and "TMBbinomial","TMBnormal", "TMBlognormal", "TMBdelta-Lognormal","TMBdelta-Gamma", "TMBnbinom1", "TMBnbinom2", and "TMBtweedie" for the corresponding models from the glmmTMB library. Binomial or TMBbinomial will be run automatically as part of the delta models if any of them are selected. @param obsdat Observer data set
#' @param obsdat  Observer data set
#' @param logdat Logbook data set
#' @param yearVar Character. The name of the year variable in \code{obsdat} and \code{logdat}. Both input files must contain the same variable name for year.
#' @param obsEffort Character. The name of the effort variable in \code{obsdat}. This variable must have the same effort units as \code{logEffort}
#' @param logEffort Character. The name of the effort variable in \code{logdat}. Optional and only used when estimating bycatch. This variable must have the same effort units as \code{obsEffort}
#' @param logUnsampledEffort Character. The name of the unsampled effort variable in in \code{logdat}. Optional and used to specify a column for effort that is not sampled, in trips with observers. This can be zero in all cases if observers sample 100% of effort in sampled trips. Only used when \code{includeObsCatch} is TRUE
#' @param includeObsCatch Logical. Set to TRUE if (1) the observed sample units can be matched to the logbook sample units and (2) you want to calculate total bycatch as the observed bycatch plus the predicted unobserved bycatch. This doesn't work with aggregated logbook effort.
#' @param matchColumn Character. If \code{includeObsCatch} is TRUE, give the name of the column that matches sample units between the observer and logbook data. Otherwise, this can be NA
#' @param factorNames Character vector. Specify which variables should be interpreted as categorical, ensuring imposes factor format on these variables. Variables not in this list will retain their original format. These variables must have identical names and factor levels in \code{obsdat} and \code{logdat}
#' @param randomEffects Character vector. Random effects that should be included in all non-delta and binomial models, as a character vector in (e.g. "Year:area" to include Year:area as a random effect). Null if none. Note that random effects will be included in all models. The code will not evaluate whether they should be included.
#' @param randomEffects2 Character vector. Random effects that should be included in the positive catch component of delta models, as a character vector in (e.g. "Year:area" to include Year:area as a random effect). Null if none. Note that random effects will be included in all models. The code will not evaluate whether they should be included.
#' @param logNum Character vector. The name of the column in \code{logdat} that gives the number of sample units (e.g., trips or sets). If the logbook data is not aggregated (i.e. each row is a sample unit) set value to NA
#' @param sampleUnit Character. What is the sample unit in \code{logdat}? e.g. sets or trips.
#' @param EstimateIndex Logical. What would you like to estimate? You may calculate either an annual abundance index, or total bycatch, or both.
#' @param EstimateBycatch Logical. What would you like to estimate? You may calculate either an annual abundance index, or total bycatch, or both. If you want total bycatch, you must provide logbook data or some other source of total effort to \code{logdat}.
#' @param complexModel Specify as stats::formula. Specify the most complex and simplest model to be considered. The code will find compare all intermediate models using information criteria. Include only fixed effects.
#' @param simpleModel Specify as stats::formula. This model includes all variables tha must be in the final bycatch estimation model
#' @param indexModel Specify as stats::formula. Use indexModel to specify which strata to keep separate in calculating abundance indices.
#' @param designMethods Character vector of methods to use for design based estimation. Current options are Ratio and Delta (for a delta-lognormal estimator).
#' @param designVars Specify strata that must be included in design based estimates, in order across which data should be pooled
#' @param designPooling TRUE if design-based estimates should be pooled for strata with missing data
#' @param poolTypes Type of pooling for each variable in designVars, as a character vector in the same order. Options are "all", "pooledVar" and (currently for year only) "adjacent"
#' @param pooledVar Variables to pool over for any variable with pooledVar in the previous line, as a character vector in the same order as designVars. Use NA for variables with other pooling methods.  This can be used to pool (for example) months into seasons when pooling is needed.
#' @param adjacentNum Number of adjacent years to include for adjacent pooling, as a numberical vector in the same order as designVars. NA for anything other than year.
#' @param minStrataUnit The smallest sample size in the strata defined by designVars that is acceptable, in sample units (e.g. trips)
#' @param baseDir Character. A directory to save output. Defaults to current working directory.
#' @param runName Characer. Give a name to the run, which will be used to set up a directory for the outputs
#' @param runDescription Character. Brief summary of the run, which will be used to set up a directory for the outputs
#' @param common Character vector. Provide a common name for the species used in bycatch and index estimation. Can be a vector of names to do multiple species at the same time.
#' @param sp  Character vector. Provide a scientific name for the species used in bycatch and index estimation. Can be a vector of names to do multiple species at the same time
#' @param obsCatch Character vector. The name of the column(s) in \code{obsdat} that contain catch. If it is a vector, order of variable names must follow the same order as names provided in \code{common} and \code{sp}
#' @param catchUnit Character vector. Give units of catch (e.g., number) to go in plot labels. Must be a vector of the same length as \code{sp}
#' @param catchType Character vector. Give type of catch (e.g., dead discards) to go in plot labels. Must be a vector of the same length as \code{sp}
#' @import ggplot2 parallel dplyr doParallel foreach utils tidyverse parallelly
#' @importFrom stats median
#' @export
#' @keywords Fitting functions
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
#' setupObj<-bycatchSetup(
#' modelTry = c("Delta-Lognormal","TMBnbinom2"),
#' obsdat = obsdatExample,
#' logdat = logdatExample,
#' yearVar = "Year",
#' obsEffort = "sampled.sets",
#' logEffort = "sets",
#' logUnsampledEffort = NULL,
#' includeObsCatch  = FALSE,
#' matchColumn = NA,
#' factorNames = c("Year","season"),
#' randomEffects = NULL,
#' randomEffects2 = NULL,
#' EstimateIndex = TRUE,
#' EstimateBycatch = TRUE,
#' logNum = NA,
#' sampleUnit = "trips",
#' complexModel = formula(y~(Year+season)^2),
#' simpleModel = formula(y~Year),
#' indexModel = formula(y~Year),
#' designMethods = c("Ratio", "Delta"),
#' designVars = c("Year","season"),
#' designPooling = TRUE,
#' poolTypes=c("adjacent","all"),
#' pooledVar=c(NA,NA),
#' adjacentNum=c(1,NA),
#' minStrataUnit = 1,
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
  randomEffects=NULL,
  randomEffects2=NULL,
  EstimateIndex=FALSE,
  EstimateBycatch,
  logNum,
  sampleUnit,
  complexModel,
  simpleModel,
  indexModel=NULL,
  designMethods = "None",
  designVars="Year",
  designPooling = FALSE,
  poolTypes=NULL,
  pooledVar=NULL,
  adjacentNum=NULL,
  minStrataUnit=1,
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
  NumCores<-parallelly::availableCores()  #Check if machine has multiple cores for parallel processing

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

  # Set up variables
  #if("Year" %in% names(obsdat) & !yearVar=="Year") obsdat<-obsdat %>% rename(oldYear=Year)
  #if("Year" %in% names(logdat) & !yearVar=="Year") logdat<-logdat %>% rename(oldYear=Year)
  obsdat<-obsdat %>% ungroup() %>%
    rename(Year=!!yearVar)
  if(EstimateBycatch) {
  logdat<-logdat %>%  ungroup() %>%
    rename(Year=!!yearVar)
  }
  requiredVarNames<-as.vector(getAllTerms(simpleModel))
  allVarNames<-as.vector(getAllTerms(complexModel))
  allVarNames<-allVarNames[grep(":",allVarNames,invert=TRUE)]
  allVarNames<-allVarNames[grep("I(*)",allVarNames,invert=TRUE)]
  allVarNames<-unique(c("Year",allVarNames))  #EAB 2/18/2025
  if(!is.null(randomEffects)) temp<-unlist(strsplit(randomEffects,":")) else temp<-NULL
  if(!is.null(randomEffects2)) temp<-c(temp,unlist(strsplit(randomEffects2,":"))) else temp<-NULL
  if(designPooling & length(pooledVar[!is.na(pooledVar)]>0)) temp2<-pooledVar[!is.na(pooledVar)] else temp2<-NULL
  allVarNames<-unique(c(allVarNames,temp,temp2,designVars))
  if(!all(allVarNames %in% names(obsdat)))
    print(paste0("Variable ", allVarNames[!allVarNames%in% names(obsdat) ], " not found in observer data"))
  if(!all(allVarNames %in% names(logdat)) & EstimateBycatch)
    print(paste0("Variable ", allVarNames[!allVarNames%in% names(logdat) ], " not found in logbook data"))
  #It's all right not to see the variable name if it is a function of another variable that is present
  if(EstimateIndex) {
   indexVarNames<-as.vector(getAllTerms(indexModel))
   if(!"Year" %in% indexVarNames) indexVarNames<-c("Year",indexVarNames)
  } else indexVarNames=NULL
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
    if(includeObsCatch & EstimateBycatch) {
      obsdat<-obsdat %>% rename(matchColumn=!!matchColumn)
      logdat<-logdat %>% rename(matchColumn=!!matchColumn,unsampledEffort=!!logUnsampledEffort)
    }
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
  #Subtract first year if numeric to improve convergence
  if(is.numeric(obsdat$Year) & "Year" %in% allVarNames) {
    startYear<-min(obsdat$Year)
    obsdat$Year<-obsdat$Year-startYear
   if(EstimateBycatch) logdat$Year<-logdat$Year-startYear
   if(EstimateIndex)  indexDat$Year<-indexDat$Year-startYear
  } else startYear<-min(as.numeric(as.character(obsdat$Year)))

  #Set up variables if looping over species
  numSp<-length(sp)
  if(numSp>1) {
    if(length(catchUnit)==1) catchUnit<-rep(catchUnit,numSp)
    if(length(catchType)==1) catchType<-rep(catchType,numSp)
  }

    #Set up directory for output
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
  yearSum<-list()
  yearSumGraph<-list()
  poolingSum<-list()
  includePool<-list()

  #Make lists to keep output, which will also be output as .pdf and .csv files for use in reports.
  dirname<-list()
  dat<-list()
  #Loop through all species and print data summary. Note that records with NA in either catch or effort are excluded automatically
  for(run in 1:numSp) {
    dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/")
    if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]])
    if(includeObsCatch & EstimateBycatch) tempvars<-c(allVarNames,"Effort","Catch","matchColumn") else
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
    #Make annual summary
    yearSum[[run]]<-MakeSummary(
                      obsdatval = dat[[run]],
                      logdatval = logdat,
                      strataVars = "Year",
                      EstimateBycatch = EstimateBycatch,
                      startYear = startYear
                    )
    if(("Ratio" %in% designMethods | "Delta" %in% designMethods) & EstimateBycatch) {
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
              paste0(dirname[[run]],common[run],catchType[run],"DataSummary.csv"), row.names = FALSE)
    if(EstimateBycatch) {
     x<-list("Unstratified ratio"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$Cat,Total.se=.data$Cse))
     if("Ratio" %in% designMethods)  x=c(x,list("Ratio"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$ratioMean,Total.se=.data$ratioSE)))
     if("Delta" %in% designMethods)  x=c(x,list("Design Delta"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$deltaMean,Total.se=.data$deltaSE)))
     yearSumGraph[[run]]<-bind_rows(x,.id="Source")     %>%
        mutate(TotalVar=.data$Total.se^2,Total.cv=.data$Total.se/.data$Total,
               Total.mean=NA,TotalLCI=.data$Total-1.96*.data$Total.se,TotalUCI=.data$Total+1.96*.data$Total.se) %>%
        mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI))
     #Calculations at level of simple model
     strataSum[[run]]<-MakeSummary(
                        obsdatval = dat[[run]],
                        logdatval = logdat,
                        strataVars = unique(c("Year",requiredVarNames)),
                        EstimateBycatch = EstimateBycatch,
                        startYear = startYear
                      )
    write.csv(strataSum[[run]],
              paste0(dirname[[run]],common[run],catchType[run],"StrataSummary.csv"), row.names = FALSE)
   }
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
                      params=list(outDir=outDir),
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
      randomEffects = randomEffects,
      randomEffects2 = randomEffects2,
      EstimateIndex = EstimateIndex,
      EstimateBycatch =EstimateBycatch,
      logNum = logNum,
      sampleUnit = sampleUnit,
      complexModel = complexModel,
      simpleModel = simpleModel,
      indexModel = indexModel,
      designMethods = designMethods,
      designVars = designVars,
      designPooling = designPooling,
      poolTypes=poolTypes,
      pooledVar=pooledVar,
      adjacentNum=adjacentNum,
      minStrataUnit = minStrataUnit,
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
      yearSumGraph = yearSumGraph,
      requiredVarNames = requiredVarNames,
      allVarNames = allVarNames,
      startYear = startYear,
      strataSum = strataSum,
      poolingSum = poolingSum,
      includePool = includePool,
      NumCores = NumCores
    )
  )

  saveRDS(output, file=paste0(outDir,"/", Sys.Date(),"_BycatchModelSpecification.rds"))
  return(output)
}

















