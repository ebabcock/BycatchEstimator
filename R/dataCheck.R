

#---------------------------------
#Data Checking
#----------------------------------

#Roxygen header
#'Make plots and tables of the data
#'
#'
#'
#' @param setupObj  An object produced by \code{bycatchSetup}.
#' @export
#' @keywords Data plotting function
#' @examples
#' \dontrun{
#' library(BycatchEstimator)
#' dataCheck(setupObj)
#' #' }
#' #-------------------------------------------------

dataCheck<-function(setupObj){
  #Setup directory naming
  outDir<-paste0(setupObj$bycatchInputs$baseDir, paste("/Output", setupObj$bycatchInputs$runName))
  if(!dir.exists(outDir)) dir.create(outDir)
  fileName<-paste0(outDir,"/", Sys.Date(),"_BycatchModelSpecification.rds")
  if(!file.exists(fileName))
    saveRDS(output, file=fileName)
  mkd<-tryCatch({
    system.file("Markdown", "printDataChecks.Rmd",
                package = "BycatchEstimator", mustWork = TRUE)
  },
  error = function(c) NULL
  )
  if(!is.null( mkd)){
    rmarkdown::render(mkd,
                      params=list(outDir=outDir),
                      output_file = "DataChecks.html",
                      output_dir=outDir,
                      quiet = TRUE)
  }


}
