
#---------------------------------
#Design-based estimators
#----------------------------------

#Roxygen header
#'Bycatch estimation using design-based estimators
#'
#' Produces estimates of bycatch using design-based ratio estimator and delta estimator, with the option of pooling across stratification variables.
#'
#'
#' @param setupObj  An object produced by \code{bycatchSetup}.
#' @param designScenario Short name, e.g. noPool, or Pool1, to distinguish outputs made with the same setupObj.
#' @param designMethods Character vector of methods to use for design based estimation. Current options are Ratio and Delta (for a delta-lognormal estimator).
#' @param designVars Specify strata that must be included in design based estimates, in order across which data should be pooled. Order of these variables determines order for which pooling will occur.
#' @param groupVar Specify variable to keep separate in summaries. Defaults to Year. Put "NA" to summarize over whole dataset
#' @param designPooling TRUE if design-based estimates should be pooled for strata with missing data
#' @param poolTypes Type of pooling for each variable in designVars, as a character vector in the same order. Options are "none", where data will not be pooled over this variable, "all" where data will be pooled over all levels of the variable, "pooledVar" where the variable named in pooledVar will be used to pool, and (currently for year only) "adjacent" to pool over adjacent years.
#' @param pooledVar Variables to pool over for any variable with pooledVar in the previous line, as a character vector in the same order as designVars. Use NA for variables with other pooling methods.  This can be used to pool (for example) months into seasons when pooling is needed.
#' @param adjacentNum Number of adjacent years to include for adjacent pooling, as a numerical vector in the same order as designVars. NA for anything other than year, in the same order as designVars.
#' @param minStrataUnit The smallest sample size in the strata defined by designVars that is acceptable, in sample units (e.g. trips); below that pooling will occur.
#' @param reportType Character. Choose type of report to be produced. Options are pdf, html (default) or both.
#' @import ggplot2 dplyr utils tidyverse
#' @importFrom stats median
#' @export
#' @keywords Fitting functions
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
#' #-------------------------------------------------
#' #Step 1. Run the datasetup function and review data inputs
#'setupObj<-bycatchSetup(
#' obsdat = obsdatExample,
#' logdat = logdatExample,
#' yearVar = "Year",
#' obsEffort = "sampled.sets",
#' logEffort = "sets",
#' obsCatch = "Catch",
#' catchUnit = "number",
#' catchType = "dead discard",
#' logNum = NA,
#' sampleUnit = "trips",
#' factorVariables = c("Year","season"),
#' numericVariables = NA,
#' EstimateBycatch = TRUE,
#' runName = "SimulatedExample",
#' runDescription = "Example with simulated data",
#' common = "Simulated species",
#' sp = "Genus species",
#' reportType = "html"
#')
#'
#' #-------------
#' #Step 2. Design-based estimators (with pooling)
#' bycatchDesign(
#' setupObj = setupObj,
#' designScenario = "noPool",
#' designMethods = c("Ratio", "Delta"),
#' designVars = c("Year","season"),
#' groupVar = "Year",
#' designPooling = TRUE,
#' poolTypes=c("adjacent","all"),
#' pooledVar=c(NA,NA),
#' adjacentNum=c(1,NA),
#' minStrataUnit = 1
#')}


bycatchDesign <- function(
  setupObj = setupObj,
  designScenario,
  designMethods = "Ratio",
  designVars = "Year",
  groupVar = "Year",
  designPooling = FALSE,
  poolTypes = NULL,
  pooledVar = NULL,
  adjacentNum = NULL,
  minStrataUnit = 1,
  reportType = "html"
){

  #unpack setup obj
  obsdat<-logdat<-yearVar<-obsEffort<-logEffort<-obsCatch<-catchUnit<-catchType<-
    logNum<-sampleUnit<-factorVariables<-numericVariables<-EstimateBycatch<-
    baseDir<-dirname<-outDir<-runName<-runDescription<-common<-sp<-NULL

  dat<-numSp<-yearSum<-allVarNames<-startYear<-strataSum<-shortName<-NULL

  for(r in 1:NROW(setupObj$bycatchInputs)) assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]]) #assign values of bycatchInputs to each element
  for(r in 1:NROW(setupObj$bycatchOutputs)) assign(names(setupObj$bycatchOutputs)[r],setupObj$bycatchOutputs[[r]])

  if(designPooling & length(pooledVar[!is.na(pooledVar)]>0)) temp2<-pooledVar[!is.na(pooledVar)] else temp2<-NULL

  #make sure only one groupVar
  groupVar<-groupVar[1]

  #check variables
  if(!all(designVars %in% c(factorVariables,"Year"))) stop(paste0("The design variables for design-based estimation must be in the list of factor variables in the setup object. Year may be a number or a factor."))


  #Set up directory for output
 #outDir<-paste0(baseDir, paste("/Output", runName))
#  if(!dir.exists(outDir)) dir.create(outDir)

  #Make R objects to store analysis
  poolingSum<-list()
  includePool<-list()
  yearSumGraph<-list()
  designyeardf <- list()
  designstratadf<-list()

  # spp loop
  for(run in 1:numSp) {
    dirname[[run]]<-gsub("Setup files","Design files",dirname[[run]])
    if(!dir.exists(dirname[[run]])) dir.create(dirname[[run]],recursive = TRUE)

  if(("Ratio" %in% designMethods | "Delta" %in% designMethods)) {
    if(designPooling) {
      temp<-getPooling(obsdatval= dat[[run]],
                       logdatval=logdat,
                       minStrataUnit=minStrataUnit,
                       designVars=designVars,
                       pooledVar=pooledVar,
                       poolTypes=poolTypes,
                       adjacentNum=adjacentNum)
      poolingSum[[run]]<-temp[[1]]
      # excluding some of columns when saving csv: needs.pooling, poolnum, stratum
      write.csv(poolingSum[[run]][,!(names(poolingSum[[run]]) %in% c("needs.pooling","stratum"))],
                paste0(dirname[[run]],shortName[run],designScenario,"Pooling.csv"),
                row.names = FALSE)
      includePool[[run]]<-temp[[2]]
    } else  {
      poolingSum[[run]]<-NULL
      includePool[[run]]<-NULL
    }

    designyeardf[[run]]<-getDesignEstimates(obsdatval = dat[[run]],
                             logdatval = logdat,
                             strataVars = groupVar,
                             designVars = designVars,
                             designPooling = designPooling,
                             minStrataUnit = minStrataUnit,
                             startYear = startYear,
                             poolingSum = poolingSum[[run]],
                             includePool= includePool[[run]]
    )
    write.csv(designyeardf[[run]],
              paste0(dirname[[run]],shortName[run],designScenario,"DesignYear.csv"),
              row.names = FALSE)

    #And design based stratification
    designstratadf[[run]]<-getDesignEstimates(obsdatval = dat[[run]],
                             logdatval = logdat,
                             strataVars = designVars,
                             designVars = designVars,
                             designPooling = designPooling,
                             minStrataUnit = minStrataUnit,
                             startYear = startYear,
                             poolingSum = poolingSum[[run]],
                             includePool= includePool[[run]]
    )
    write.csv(designstratadf[[run]],
              paste0(dirname[[run]],shortName[run],designScenario,"DesignStrata.csv"), row.names = FALSE)
  }
  if(all(is.na(groupVar))) {
    yearSum[[run]]$Year="All"
    designyeardf[[run]]$Year="All"
  }
  if(groupVar!="Year" & !is.na(groupVar)) {
    x<-NULL
    if("Ratio" %in% designMethods)
      x=c(x,list("Ratio"=dplyr::select(designyeardf[[run]],!!groupVar,Total=.data$ratioMean,Total.se=.data$ratioSE)))
    if("Delta" %in% designMethods)
      x=c(x,list("Design Delta"=dplyr::select(designyeardf[[run]],!!groupVar,Total=.data$deltaMean,Total.se=.data$deltaSE)))
  }
    if(groupVar=="Year") {
      x<-list("Unstratified ratio"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$Cat,Total.se=.data$Cse))
      if("Ratio" %in% designMethods)
        x=c(x,list("Ratio"=dplyr::select(designyeardf[[run]],Year=.data$Year,Total=.data$ratioMean,Total.se=.data$ratioSE)))
      if("Delta" %in% designMethods)
        x=c(x,list("Design Delta"=dplyr::select(designyeardf[[run]],Year=.data$Year,Total=.data$deltaMean,Total.se=.data$deltaSE)))
    }

    yearSumGraph[[run]]<-bind_rows(x,.id="Source")     %>%
      mutate(TotalVar=.data$Total.se^2,Total.cv=.data$Total.se/.data$Total,
             Total.mean=NA,TotalLCI=.data$Total-1.96*.data$Total.se,TotalUCI=.data$Total+1.96*.data$Total.se) %>%
      mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI))

  } #close loop for each sp


  #Create output list
  output<-list(

    designInputs = list(
      designScenario=designScenario,
      designMethods = designMethods,
      designVars = designVars,
      groupVar = groupVar,
      designPooling = designPooling,
      poolTypes=poolTypes,
      pooledVar=pooledVar,
      adjacentNum=adjacentNum,
      minStrataUnit = minStrataUnit,
      baseDir = baseDir
    ),

    designOutputs = list(
      yearSum = yearSum,
      yearSumGraph = yearSumGraph,
      strataSum = strataSum,
      poolingSum = poolingSum,
      includePool = includePool,
      designyeardf = designyeardf,
      designstratadf = designstratadf
    )

  )

  saveRDS(output, file=paste0(outDir,"/", Sys.Date(),"_BycatchDesign",designScenario,".rds"))


  #Create report
  for(run in 1:numSp) {
    if(reportType == "html" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "printBycatchDesign.Rmd",
                    package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )

      if(!is.null( mkd)){

        rmarkdown::render(mkd,
                          params=list(outDir=outDir, run = run, designScenario=designScenario),
                          output_format = "html_document",
                          output_file = paste0(shortName[run],designScenario ,"DesignResults.html"),
                          output_dir=paste0(outDir,"/",shortName[run],"/"),
                          quiet = TRUE)
      }
    }

    if(reportType == "pdf" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "printBycatchDesign.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )

      if(!is.null( mkd)){

        tryCatch({
        rmarkdown::render(mkd,
                          params=list(outDir=outDir, run = run),
                          output_format = "pdf_document",
                          output_file = paste0(shortName[run],designScenario , "DesignResults.pdf"),
                          output_dir=paste0(outDir,"/",shortName[run],"/"),
                          quiet = TRUE)

        },
        error = function(e){
          message("PDF rendering failed, reverting to html.")
          rmarkdown::render(mkd,
                            params=list(outDir=outDir, run = run),
                            output_format = "html_document",
                            output_file = paste0(shortName[run],designScenario ,"DesignResults.html"),
                            output_dir=paste0(outDir,"/",shortName[run],"/"),
                            quiet = TRUE)
          })
        }


      }

   }

      #Clean up: delete the figures/ directory after rendering
      if(!is.null(mkd)) {
        fig_dir <- file.path(dirname(mkd), "figures")
        if (dir.exists(fig_dir)) {
          unlink(fig_dir, recursive = TRUE)
        }
      }


  #return(output)

} #close main function



