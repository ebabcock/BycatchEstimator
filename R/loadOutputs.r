
#---------------------------------
#Load R objects for further analysis
#----------------------------------

#Roxygen header
#' loadOutputs
#'
#' Reads in all the R objects created during runs of the bycatchEstimator for further analysis
#'
#' @param baseDir The base directory for the runs, same as bycatchSetup.
#' @param runName The run name, same as bycatchSetup.
#' @param runDate  The date when the model was run. Defaults to current date, but can be set to read in models previously run.
#' @param designScenarios Character vector of designScenario values from original run. NULL to read in no design-based results.
#' @param modelScenarios Character vector of modelScenario values from original run.
#' @export
#' @returns Returns a list with the setupObj from the specified run, a list called designObjList
#' which contains the design-based model inputs and outputs for each designScenario, modelObjList,
#'  which is the same for the models in modelScenarios, a data frame called allYearEstimates which
#'  is the annual estimates across all design-based and model-based scenarios in a format suitable
#'   for ggplot, allGroupEstimates, which is the same by the groups defined in groupVar
#'   and predictionGroups, as well as the results by groups for models and design-based estimators
#'   separately (allModResultsGroup, allDesignResultsGroup).
#' @keywords reload outputs
loadOutputs<-function(baseDir = getwd(),
                      runName,
                      runDate =  Sys.Date(),
                      designScenarios = NULL,
                      modelScenarios = NULL) {
  #Check that setup file exists and read in.
  tempRunName<-abbreviate(runName,minlength=10)
  outDir<-paste0(baseDir, paste0("/Output", tempRunName))
  if(!dir.exists(outDir)) stop(paste("Directory",outDir,"not found."))
  setupFile<-paste0(runDate,"_BycatchSetupSpecification.rds")
  if(!file.exists(paste0(outDir,"/",setupFile))) stop(paste("Setup file",setupFile ,"not found in",outDir,"."))
  setupObj<-readRDS(file=paste0(outDir,"/",runDate,"_BycatchSetupSpecification.rds"))
  numSp<-setupObj$bycatchOutputs$numSp
  if(all(is.na(designScenarios))) designScenarios<-NULL
  #If doing design based, check that design file exists and read in.
  if(!is.null(designScenarios)) {
    designObjList<-list()
    allDesignResults<-list()
    allDesignResultsGroup<-list()
    for(i in 1:length(designScenarios)) {
      designFile<-paste0(runDate,"_BycatchDesign",designScenarios[i],".rds")
      if(!file.exists(paste0(outDir,"/",designFile))) stop(paste("Design file",designFile ,"not found in",outDir,"."))
      designObjList[[i]]<-readRDS(file=paste0(outDir,"/",designFile))
      if(length(designObjList[[i]]$designOutputs$yearSumGraph)>0) {
        names(designObjList[[i]]$designOutputs$yearSumGraph)<-
          paste(1:numSp,setupObj$bycatchInputs$common,setupObj$bycatchInputs$sp,setupObj$bycatchInputs$catchType,sep=";")
        allDesignResults[[i]]<-bind_rows(designObjList[[i]]$designOutputs$yearSumGraph,
                                         .id="Common")%>%
          separate_wider_delim(Common,delim=";",names = c("spNum","Common","Species","CatchType"))%>%
          mutate(Valid=1)
        names(allDesignResults)[i]<-designScenarios[i]
        if("Year" %in% names(allDesignResults[[i]]))
          if(!is.numeric(allDesignResults[[i]]$Year))
            allDesignResults[[i]]$Year<-as.numeric(as.character(allDesignResults[[i]]$Year))
      }
      if(length(designObjList[[i]]$designOutputs$groupSumGraph)>0) {
        names(designObjList[[i]]$designOutputs$groupSumGraph)<-
          paste(1:numSp,setupObj$bycatchInputs$common,setupObj$bycatchInputs$sp,setupObj$bycatchInputs$catchType,sep=";")
        allDesignResultsGroup[[i]]<-bind_rows(designObjList[[i]]$designOutputs$groupSumGraph,
                                               .id="Common")%>%
          separate_wider_delim(Common,delim=";",names = c("spNum","Common","Species","CatchType"))%>%
          mutate(Valid=1)
        names(allDesignResultsGroup)[i]<-designScenarios[i]
        if("Year" %in% names(allDesignResultsGroup[[i]]))
          if(!is.numeric(allDesignResultsGroup[[i]]$Year))
            allDesignResultsGroup[[i]]$Year<-as.numeric(as.character(allDesignResultsGroup[[i]]$Year))
      }
    }
    if(length(allDesignResults) > 0)
      allDesignResults<-bind_rows(allDesignResults,.id="Scenario")
    if(length(allDesignResultsGroup)>0)
      allDesignResultsGroup<-bind_rows(allDesignResultsGroup,.id="Scenario")
    if("Year" %in% names(allDesignResults)) {
      allDesignResults<-mutate(allDesignResults,Year=as.numeric(as.character(Year)))
      allDesignResultsGroup<-mutate(allDesignResultsGroup,Year=as.numeric(as.character(Year)))
    }
  } else {
    allDesignResults<-NULL
    designObjList<-NULL
  }
  if(all(is.na(modelScenarios))) modelScenarios<-NULL
  if(!is.null(modelScenarios)) {
    modelObjList<-list()
    allModResults<-list()
    allModResultsGroup<-list()
    for(i in 1:length(modelScenarios)) {
      modelFile<-paste0(runDate,"_BycatchFit",modelScenarios[i],".rds")
      if(!file.exists(paste0(outDir,"/",modelFile))) stop(paste("Model file",modelFile ,"not found in",outDir,"."))
      modelObjList[[i]]<-readRDS(file=paste0(outDir,"/",modelFile))
      names(modelObjList[[i]]$modelOutputs$allmods)<-paste(1:numSp,setupObj$bycatchInputs$common,setupObj$bycatchInputs$sp,setupObj$bycatchInputs$catchType,sep=";")
      allModResults[[i]]<-bind_rows(modelObjList[[i]]$modelOutputs$allmods,
                                    .id="Common")%>%
        separate_wider_delim(Common,delim=";",names = c("spNum","Common","Species","CatchType"))
      names(modelObjList[[i]]$modelOutputs$allmodsStrata)<-paste(1:numSp,setupObj$bycatchInputs$common,setupObj$bycatchInputs$sp,setupObj$bycatchInputs$catchType,sep=";")
      allModResultsGroup[[i]]<-bind_rows(modelObjList[[i]]$modelOutputs$allmodsStrata,
                                          .id="Common")%>%
        separate_wider_delim(Common,delim=";",names = c("spNum","Common","Species","CatchType"))
    }
    names(modelObjList)<-modelScenarios
    names(allModResults)<-modelScenarios
    allModResults<-bind_rows(allModResults,.id="Scenario")
    allModResultsGroup<-bind_rows(allModResultsGroup,.id="Scenario")
    if("Year" %in% names(allModResults)) {
      allModResults<-mutate(allModResults,Year=as.numeric(as.character(Year)))
      allModResultsGroup<-mutate(allModResultsGroup,Year=as.numeric(as.character(Year)))
    }
  }  else {
    allModResults<-NULL
    allModResultsGroup<-NULL
    modelObjList<-NULL
  }
  allYearEstimates<-bind_rows(allModResults,allDesignResults) %>%
    mutate(Run=runName)
  if(nrow(allYearEstimates)>0) allYearEstimates<-filter(allYearEstimates,!Source=="Unstratified ratio")
  if(length(setdiff(names(allDesignResultsGroup),names(allModResultsGroup)))==0) {
    allGroupEstimates <-
      bind_rows(allDesignResultsGroup,allModResultsGroup) %>%
      mutate(Run=runName)
  } else
    allGroupEstimates<-NULL
    list(setupObj=setupObj,
       designObjList=designObjList,
       modelobjList=modelObjList,
       allYearEstimates=allYearEstimates,
       allGroupEstimates=allGroupEstimates,
       allModResultsGroup=allModResultsGroup,
       allDesignResultsGroup=allDesignResultsGroup,
       runName=runName,
       baseDir=baseDir,
       runDate=runDate
  )

}
