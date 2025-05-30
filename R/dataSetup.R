
#---------------------------------
#Data setup (new)
#----------------------------------

#Roxygen header
#'Bycatch estimation data setup (new)
#'
#'Sets global conditions, makes a preliminary data summary and data checks (make plots and tables of the data).
#'
#'
#' @param obsdat  Observer data set
#' @param logdat Logbook data set
#' @param yearVar Character. The column name of the year variable in \code{obsdat} and \code{logdat}. Both input files must contain the same variable name for year.
#' @param obsEffort Character. The column name of the effort variable in \code{obsdat}. This variable must have the same effort units as \code{logEffort}
#' @param logEffort Character. The column name of the effort variable in \code{logdat}. Optional and only used when estimating bycatch. This variable must have the same effort units as \code{obsEffort}
#' @param obsCatch Character vector. The name of the column(s) in \code{obsdat} that contain catch. If it is a vector, order of variable names must follow the same order as names provided in \code{common} and \code{sp}
#' @param catchUnit Character vector. Give units of catch (e.g., number) to go in plot labels. Must be a vector of the same length as \code{sp}
#' @param catchType Character vector. Give type of catch (e.g., dead discards) to go in plot labels. Must be a vector of the same length as \code{sp}
#' @param logNum Character vector. The name of the column in \code{logdat} that gives the number of sample units (e.g., trips or sets). If the logbook data is not aggregated (i.e. each row is a sample unit) set value to NA
#' @param sampleUnit Character. What is the sample unit in \code{logdat}? e.g. sets or trips.
#' @param factorVariables Character vector. Specify which variables should be interpreted as categorical, ensuring imposes factor format on these variables. These variables must have identical names and factor levels in \code{obsdat} and \code{logdat}
#' @param numericVariables Character vector. Specify which variables should be interpreted as numeric. These variables must have identical names in \code{obsdat} and \code{logdat}. If there are no numeric variables, set numericVariables=NA.
#' @param includeObsCatch Logical. Set to TRUE if (1) the observed sample units can be matched to the logbook sample units and (2) you want to calculate total bycatch as the observed bycatch plus the predicted unobserved bycatch. This doesn't work with aggregated logbook effort.
#' @param logUnsampledEffort Character. The name of the unsampled effort variable in in \code{logdat}. Optional and used to specify a column for effort that is not sampled, in trips with observers. This can be zero in all cases if observers sample 100% of effort in sampled trips. Only used when \code{includeObsCatch} is TRUE
#' @param matchColumn Character. If \code{includeObsCatch} is TRUE, give the name of the column that matches sample units between the observer and logbook data. Otherwise, this can be NA
#' @param EstimateIndex Logical. What would you like to estimate? You may calculate either an annual abundance index, or total bycatch, or both.
#' @param EstimateBycatch Logical. What would you like to estimate? You may calculate either an annual abundance index, or total bycatch, or both. If you want total bycatch, you must provide logbook data or some other source of total effort to \code{logdat}.
#' @param baseDir Character. A directory to save output. Defaults to current working directory.
#' @param runName Characer. Give a name to the run, which will be used to set up a directory for the outputs
#' @param runDescription Character. Brief summary of the run, which will be used to set up a directory for the outputs
#' @param common Character vector. Provide a common name for the species used in bycatch and index estimation. Can be a vector of names to do multiple species at the same time.
#' @param sp  Character vector. Provide a scientific name for the species used in bycatch and index estimation. Can be a vector of names to do multiple species at the same time
#' @param reportType Character. Choose type of data checks report to be produced. Options are pdf, html (default) or both.
#' @import ggplot2 parallel dplyr doParallel foreach utils tidyverse parallelly
#' @importFrom stats median
#' @export
#' @keywords Fitting functions
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
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
#' EstimateBycatch = TRUE,
#' baseDir = getwd(),
#' runName = "SimulatedExample",
#' runDescription = "Example with simulated data",
#' common = "Simulated species",
#' sp = "Genus species",
#' reportType = "html"
#')}

bycatchSetup_new <- function(
    obsdat,
    logdat,
    yearVar,
    obsEffort,
    logEffort,
    obsCatch,
    catchUnit,
    catchType,
    logNum,
    sampleUnit,
    factorVariables,
    numericVariables,
    includeObsCatch  = FALSE,
    logUnsampledEffort = NULL,
    matchColumn = NA,
    EstimateIndex=FALSE,
    EstimateBycatch=TRUE,
    baseDir = getwd(),
    runName,
    runDescription,
    common,
    sp,
    reportType = "html"
){

  SampleUnits<-Year<-drop_na<-Catch<-Effort<-cpue<-pres<-Pos<-OUnit<-OEff<-Eff<-Units<-outDir<-NULL
    # drop_na: drops rows with missing values

  #Set global conditions
  theme_set(theme_bw()) #ggplot theme

  defaultW <- 0 # allow code to run with warnings showing as they happen
  options(warn=defaultW)

  # Set up variables
  obsdat<-obsdat %>% ungroup() %>%
    rename(Year=!!yearVar) # !! unquotes variable and extracts value, not name of the variable
  if(EstimateBycatch) {
    logdat<-logdat %>%  ungroup() %>%
      rename(Year=!!yearVar)
  }

  # set up factor variables and numeric variables
  if(unique(!is.na(numericVariables))){
    allVarNames <- c(factorVariables,numericVariables)}
  else{
    allVarNames <- c(factorVariables)}
  allVarNames <- unique(c("Year",allVarNames))

  if(!all(allVarNames %in% names(obsdat)))
    print(paste0("Variable ", allVarNames[!allVarNames%in% names(obsdat) ], " not found in observer data"))

  if(!all(allVarNames %in% names(logdat)) & EstimateBycatch)
    print(paste0("Variable ", allVarNames[!allVarNames%in% names(logdat) ], " not found in logbook data"))
  #It's all right not to see the variable name if it is a function of another variable that is present


  #Set up data frames
  obsdat<-obsdat %>%
    rename(Effort=!!obsEffort) %>%
    mutate(across(all_of(factorVariables),factor))

  # add here warning message about NAs in observer data
  if(any(is.na(obsdat[,c(allVarNames,Effort)]))){
    warning("There are NAs in observer data in either factor variables, numeric variables or effort. Check what to do with these before continuing.")
  }

  if(EstimateBycatch) {
    if(is.na(logNum))   {
      logdat<-mutate(logdat,SampleUnits=1)
      logNum<-"SampleUnits"
    }
    logdat<-logdat %>%
      rename(Effort=!!logEffort,SampleUnits=!!logNum) %>%
      mutate(across(all_of(factorVariables),factor))

    if(logEffort==sampleUnit) logdat<-mutate(logdat,Effort=SampleUnits)

    # add here warning message about NAs in logbook data
    if(any(is.na(logdat[,c(allVarNames,Effort)]))){
      warning("There are NAs in logbok data in either factor variables, numeric variables or effort. Check what to do with these before continuing.")
    }

    if(includeObsCatch & EstimateBycatch) {
      obsdat<-obsdat %>% rename(matchColumn=!!matchColumn)
      logdat<-logdat %>% rename(matchColumn=!!matchColumn,unsampledEffort=!!logUnsampledEffort)

      missing_trips <- setdiff(obsdat$matchColumn,logdat$matchColumn)
      if(length(missing_trips)>0){
        stop(paste("The following sample units are missing in the logbook data: "),
                paste(missing_trips,collapse = ", "))}
        }

    # add here warning message about missing factor levels
    # add here warning message about ranges of numeric variables not matching between observer and logbook data

      }


  # define startYear
  if(is.numeric(obsdat$Year) & "Year" %in% allVarNames) {
    startYear<-min(obsdat$Year)
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
  strataSum<-list()
  yearSum<-list()


  #Make lists to keep output, which will also be output as .pdf and .csv files for use in reports.
  dirname<-list()
  dat<-list()

  #Loop through all species and print data summary and data checks. Note that records with NA in either catch or effort are excluded automatically
  for(run in 1:numSp) {
    dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/","bycatchSetup files/")
    if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]],recursive = TRUE)

    if(includeObsCatch & EstimateBycatch) tempvars<-c(allVarNames,"Effort","Catch","matchColumn") else
      tempvars<-c(allVarNames,"Effort","Catch")

    dat[[run]]<-obsdat %>%
      rename(Catch=!!obsCatch[run])%>%
      dplyr::select_at(all_of(tempvars)) %>%
      drop_na()   %>%
      mutate(cpue=Catch/Effort,
             log.cpue=log(Catch/Effort),
             pres=ifelse(cpue>0,1,0))

    if(dim(dat[[run]])[1]<dim(obsdat)[1]) print(paste0("Removed ",dim(obsdat)[1]-dim(dat[[run]])[1]," rows with NA values for ",common[run]))

    #if(EstimateBycatch) { # move this if statement, because can do annual summary from observer data without EstimateBycatch
    #Make annual summary
    yearSum[[run]]<-MakeSummary(
      obsdatval = dat[[run]],
      logdatval = logdat,
      strataVars = "Year", #strataVars argument within MakeSummary function, not to be changed
      EstimateBycatch = EstimateBycatch,
      startYear = startYear
    )
    write.csv(yearSum[[run]],
              paste0(dirname[[run]],common[run]," ",catchType[run]," DataSummary.csv"), row.names = FALSE)


      #Calculations at level of simple model
      strataSum[[run]]<-MakeSummary(
        obsdatval = dat[[run]],
        logdatval = logdat,
        strataVars = unique(c("Year",factorVariables)),
        EstimateBycatch = EstimateBycatch,
        startYear = startYear
      )
      write.csv(strataSum[[run]],
                paste0(dirname[[run]],common[run]," ",catchType[run]," StrataSummary.csv"), row.names = FALSE)
    #}

  } #close loop for each spp

  #Create output list
  output<-list(

    #Inputs to bycatchDesign and bycatchFit
    bycatchInputs = list(
      obsdat = obsdat,
      logdat = logdat,
      yearVar = yearVar,
      obsEffort = obsEffort,
      logEffort = logEffort,
      obsCatch = obsCatch,
      catchUnit = catchUnit,
      catchType = catchType,
      logNum = logNum,
      sampleUnit = sampleUnit,
      factorVariables = factorVariables,
      numericVariables = numericVariables,
      logUnsampledEffort = logUnsampledEffort,
      includeObsCatch  = includeObsCatch,
      matchColumn = matchColumn,
      EstimateIndex = EstimateIndex,
      EstimateBycatch =EstimateBycatch,
      baseDir = baseDir,
      runName = runName,
      runDescription = runDescription,
      common = common,
      sp = sp
    ),

    #Outputs from data processing
    bycatchOutputs = list(
      dat = dat,
      numSp = numSp,
      yearSum = yearSum,
      allVarNames = allVarNames,
      startYear = startYear,
      strataSum = strataSum
    )
  )

  saveRDS(output, file=paste0(outDir,"/", Sys.Date(),"_BycatchSetupSpecification.rds"))

  #Create report
  for(run in 1:numSp) {
    if(reportType == "html" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "PrintBycatchSetup.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
        error = function(c) NULL
        )

      if(!is.null( mkd)){

      rmarkdown::render(mkd,
                      params=list(outDir=outDir, run = run),
                      output_format = "html_document",
                      output_file = paste0(common[run], " ",catchType[run], " Data checks.html"),
                      output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                      quiet = TRUE)

      }
    }

    if(reportType == "pdf" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "PrintBycatchSetup.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )

      if(!is.null( mkd)){

      tryCatch({
      rmarkdown::render(mkd,
                      params=list(outDir=outDir, run = run),
                      output_format = "pdf_document",
                      output_file = paste0(common[run], " ",catchType[run], " Data checks.pdf"),
                      output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                      quiet = TRUE)
      },
      error = function(e){
        message("PDF rendering failed, reverting to html.")
        rmarkdown::render(mkd,
                          params=list(outDir=outDir, run = run),
                          output_format = "html_document",
                          output_file = paste0(common[run], " ",catchType[run], " Data checks.html"),
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

}



