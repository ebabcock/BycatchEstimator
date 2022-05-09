

#----------------------------------------------------------------------------------------------
# LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022
#----------------------------------------------------------------------------------------------
library(BycatchEstimator)
library(MuMIn)

#-------------------------------------------------
#Step 1. Run the setup file and review data inputs

setupObj<-bycatchSetup(
  modelTry = c("Lognormal","Delta-Lognormal","Delta-Gamma","TMBnbinom1","TMBnbinom2","TMBtweedie"),
  obsdat = LLSIM_BUM_Example_observer,
  logdat = LLSIM_BUM_Example_logbook,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  logUnsampledEffort = "unsampledEffort",
  includeObsCatch  = TRUE,
  matchColumn = "trip",
  factorNames = c("Year","fleet","area","season"),
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  logNum = NA,
  sampleUnit = "trips",
  complexModel = formula(y~Year+fleet+hbf+area+season+Year:area),
  simpleModel = formula(y~Year+fleet+area),
  indexModel = formula(y~Year+area),
  baseDir = getwd(),
  runName = "LLSIMBUMtrip2022Aprilobs05mc",
  runDescription = "LLSIm BUM by trip, with 5% observer coverage including observed catch in totals April 17 2022",
  common = c("Swordfish","Blue marlin")[2],
  sp = c("Xiphias gladius","Makaira nigricans")[2],
  obsCatch = c("SWO","BUM")[2],
  catchUnit = "number",
  catchType = "catch"
)


#----------------------------------------------------------------------------
#Optionally, instead of running setup, you can read-in a previous model setup

setupObj<-readRDS(file=paste0(getwd(), paste("/Output", "LLSIMBUMtrip2022Aprilobs05mc"),"/", "2022-05-02","_BycatchModelSpecification.rds"))



#-------------
#Step 2. Model Fitting

bycatchFit(
  setupObj = setupObj,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  DredgeCrossValidation = FALSE,
  ResidualTest = FALSE,
  CIval = 0.05,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 1000,
  baseDir = getwd(),
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL
)






