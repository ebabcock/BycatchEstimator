
#---------------------------------
#Design-based estimators (new)
#----------------------------------

#Roxygen header
#'Bycatch estimation using design-based estimators
#'
#' Produces estimates of bycatch using design-based ratio estimator and delta estimator, with the option of pooling across stratification variables.
#'
#'
#' @param setupObj  An object produced by \code{bycatchSetup_new}.
#' @param designMethods Character vector of methods to use for design based estimation. Current options are Ratio and Delta (for a delta-lognormal estimator).
#' @param designVars Specify strata that must be included in design based estimates, in order across which data should be pooled. Order of these variables determines order for which pooling will occur.
#' @param designPooling TRUE if design-based estimates should be pooled for strata with missing data
#' @param poolTypes Type of pooling for each variable in designVars, as a character vector in the same order. Options are "all", "pooledVar" and (currently for year only) "adjacent"
#' @param pooledVar Variables to pool over for any variable with pooledVar in the previous line, as a character vector in the same order as designVars. Use NA for variables with other pooling methods.  This can be used to pool (for example) months into seasons when pooling is needed.
#' @param adjacentNum Number of adjacent years to include for adjacent pooling, as a numerical vector in the same order as designVars. NA for anything other than year, in the same order as designVars.
#' @param minStrataUnit The smallest sample size in the strata defined by designVars that is acceptable, in sample units (e.g. trips); below that pooling will occur.
#' @param baseDir Character. A directory to save output. Defaults to current working directory.
#' @param reportType Character. Choose type of report with results to be produced. Options are pdf, html (default) or both.
#' @import ggplot2 dplyr utils tidyverse
#' @importFrom stats median
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
#' numericVariables = NA,
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
#' #Step 2. Design-based estimators (with pooling)
#'bycatchDesign(
#' setupObj = setupObj,
#' designMethods = c("Ratio", "Delta"),
#' designVars = c("Year","season"),
#' designPooling = TRUE,
#' poolTypes=c("adjacent","all"),
#' pooledVar=c(NA,NA),
#' adjacentNum=c(1,NA),
#' minStrataUnit = 1
#')}


bycatchDesign_new <- function(
  setupObj = setupObj,
  designMethods = "Ratio",
  designVars = "Year",
  designPooling = FALSE,
  poolTypes = NULL,
  pooledVar = NULL,
  adjacentNum = NULL,
  minStrataUnit = 1,
  baseDir = getwd(),
  reportType = "html"
){

  #unpack setup obj
  obsdat<-logdat<-yearVar<-obsEffort<-logEffort<-obsCatch<-catchUnit<-catchType<-
    logNum<-sampleUnit<-factorVariables<-numericVariables<-
    logUnsampledEffort<-includeObsCatch<-matchColumn<-
    EstimateIndex<-EstimateBycatch<-
    baseDir<-runName<-runDescription<-common<-sp<-NULL

  dat<-numSp<-yearSum<-allVarNames<-startYear<-strataSum<-NULL

  for(r in 1:NROW(setupObj$bycatchInputs)) assign(names(setupObj$bycatchInputs)[r], setupObj$bycatchInputs[[r]]) #assign values of bycatchInputs to each element
  for(r in 1:NROW(setupObj$bycatchOutputs)) assign(names(setupObj$bycatchOutputs)[r],setupObj$bycatchOutputs[[r]])

  if(designPooling & length(pooledVar[!is.na(pooledVar)]>0)) temp2<-pooledVar[!is.na(pooledVar)] else temp2<-NULL

  #Set up directory for output
  outDir<-paste0(baseDir, paste("/Output", runName))
  if(!dir.exists(outDir)) dir.create(outDir)

  #Make R objects to store analysis
  poolingSum<-list()
  includePool<-list()
  yearSumGraph<-list()
  dirname <- list()
  designyeardf <- list()
  designstratadf<-list()

  # spp loop
  for(run in 1:numSp) {
    dirname[[run]]<-paste0(outDir,"/",common[run]," ",catchType[run],"/","bycatchDesign files/")
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
      # exclude some of columns when saving csv: needs.pooling, poolnum, stratum
      write.csv(poolingSum[[run]][,!(names(poolingSum[[run]]) %in% c("needs.pooling","stratum"))],
                paste0(dirname[[run]],common[run],catchType[run],"Pooling.csv"), row.names = FALSE)
      includePool[[run]]<-temp[[2]]
    } else  {
      poolingSum[[run]]<-NULL
      includePool[[run]]<-NULL
    }

    designyeardf[[run]]<-getDesignEstimates(obsdatval = dat[[run]],
                             logdatval = logdat,
                             strataVars = "Year",
                             designVars = designVars,
                             designPooling = designPooling,
                             minStrataUnit = minStrataUnit,
                             startYear = startYear,
                             poolingSum = poolingSum[[run]],
                             includePool= includePool[[run]]
    )
    write.csv(designyeardf[[run]],
              paste0(dirname[[run]],common[run],catchType[run],"DesignYear.csv"), row.names = FALSE)

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
              paste0(dirname[[run]],common[run],catchType[run],"DesignStrata.csv"), row.names = FALSE)
  }


    x<-list("Unstratified ratio"=dplyr::select(yearSum[[run]],Year=.data$Year,Total=.data$Cat,Total.se=.data$Cse))
    if("Ratio" %in% designMethods)
      x=c(x,list("Ratio"=dplyr::select(designyeardf[[run]],Year=.data$Year,Total=.data$ratioMean,Total.se=.data$ratioSE)))
    if("Delta" %in% designMethods)
      x=c(x,list("Design Delta"=dplyr::select(designyeardf[[run]],Year=.data$Year,Total=.data$deltaMean,Total.se=.data$deltaSE)))

    yearSumGraph[[run]]<-bind_rows(x,.id="Source")     %>%
      mutate(TotalVar=.data$Total.se^2,Total.cv=.data$Total.se/.data$Total,
             Total.mean=NA,TotalLCI=.data$Total-1.96*.data$Total.se,TotalUCI=.data$Total+1.96*.data$Total.se) %>%
      mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI)) #this is only being used to plot all bycatch estimation methods together in one plot

  } #close loop for each sp


  #Create output list
  output<-list(

    designInputs = list(
      designMethods = designMethods,
      designVars = designVars,
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

  saveRDS(output, file=paste0(outDir,"/", Sys.Date(),"_BycatchDesignSpecification.rds"))


  #Create report
  for(run in 1:numSp) {
    if(reportType == "html" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "PrintBycatchDesign.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )

      if(!is.null( mkd)){

        rmarkdown::render(mkd,
                          params=list(outDir=outDir, run = run),
                          output_format = "html_document",
                          output_file = paste0(common[run], " ",catchType[run], " Design-based estimation results.html"),
                          output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                          quiet = TRUE)
      }
    }

    if(reportType == "pdf" || reportType == "both"){

      mkd<-tryCatch({
        system.file("Markdown", "PrintBycatchDesign.Rmd", package = "BycatchEstimator", mustWork = TRUE)
      },
      error = function(c) NULL
      )

      if(!is.null( mkd)){

        tryCatch({
        rmarkdown::render(mkd,
                          params=list(outDir=outDir, run = run),
                          output_format = "pdf_document",
                          output_file = paste0(common[run], " ",catchType[run], " Design-based estimation results.pdf"),
                          output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                          quiet = TRUE)

        },
        error = function(e){
          message("PDF rendering failed, reverting to html.")
          rmarkdown::render(mkd,
                            params=list(outDir=outDir, run = run),
                            output_format = "html_document",
                            output_file = paste0(common[run], " ",catchType[run], " Design-based estimation results.html"),
                            output_dir=paste0(outDir,"/",common[run]," ",catchType[run],"/"),
                            quiet = TRUE)
        })
      }

      #Clean up: delete the figures/ directory after rendering
      fig_dir <- file.path(dirname(mkd), "figures")
      if (dir.exists(fig_dir)) {
        unlink(fig_dir, recursive = TRUE)
      }

    }
  }

  return(output)

} #close main function



