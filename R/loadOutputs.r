
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
#' @param loadDesign  TRUE to read in design-based estimator results.
#' @param designScenario Value of designScenario from original run.
#' @param loadModel TRUE to read in model-based estimator results.
#' @param modelScenario Value of designScenario from original run.
#' @export
#' @keywords reload outputs
loadOutputs<-function(baseDir = getwd(),
                      runName,
                      runDate =  Sys.Date(),
                      loadDesign = TRUE,
                      designScenario = NULL,
                      loadModel = TRUE,
                      modelScenario = NULL) {
  #Check that setup file exists and read in.
  tempRunName<-abbreviate(runName,minlength=10)
  outDir<-paste0(baseDir, paste0("/Output", tempRunName))
  if(!dir.exists(outDir)) stop(paste("Directory",outDir,"not found."))
  setupFile<-paste0(runDate,"_BycatchSetupSpecification.rds")
  if(!file.exists(paste0(outDir,"/",setupFile))) stop(paste("Setup file",setupFile ,"not found in",outDir,"."))
  setupObj<-readRDS(file=paste0(outDir,"/",runDate(),"_BycatchSetupSpecification.rds"))
  list2env(setupObj$bycatchInputs, envir = .GlobalEnv)
  list2env(setupObj$bycatchOutputs, envir = .GlobalEnv)
  #If doing design based, check that design file exists and read in.
  if(loadDesign) {
    designFile<-paste0(runDate,"_BycatchDesign",designScenario,".rds")
    if(!file.exists(paste0(outDir,"/",designFile))) stop(paste("Design file",designFile ,"not found in",outDir,"."))
    designObj<-readRDS(file=paste0(outDir,"/",designFile))
    list2env(designObj$designInputs, envir = .GlobalEnv)
    list2env(designObj$designOutputs, envir = .GlobalEnv)
  }
  if(loadModel) {
    modelFile<-paste0(runDate,"_BycatchFit",modelScenario,".rds")
    if(!file.exists(paste0(outDir,"/",modelFile))) stop(paste("Model file",modelFile ,"not found in",outDir,"."))
    modelObj<-readRDS(file=paste0(outDir,"/",modelFile))
    list2env(modelObj$modelInputs, envir = .GlobalEnv)
    list2env(modelObj$modelOutputs, envir = .GlobalEnv)
  }
}

