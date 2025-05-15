devtools::load_all()

obsdatExample <- obsdatExample
logdatExample <- logdatExample

working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

setupObj_new <- bycatchSetup_new(
  obsdat = obsdatExample,
  logdat = logdatExample,
  yearVar = "Year",
  obsEffort = "sampled.sets",
  logEffort = "sets",
  obsCatch = "Catch",
  catchUnit = "number",
  catchType = "discard",
  logNum = NA,
  sampleUnit = "trips",
  factorVariables = c("Year","season"),
  numericVariables = NA,
  baseDir = working.dir,
  runName = "SimulatedExample_new3functions",
  runDescription = "SimulatedExample_new3functions",
  common = "Simulated species",
  sp = "Genus species"
)

design_est <- bycatchDesign_new(
  setupObj = setupObj_new,
  designMethods = c("Ratio", "Delta"),
  designVars = c("Year","season"),
  designPooling = TRUE,
  poolTypes=c("adjacent","all"),
  pooledVar=c(NA,NA),
  adjacentNum=c(1,NA),
  minStrataUnit = 1
  #baseDir = working.dir
)

model_est <- bycatchFit_new(
  setupObj = setupObj_new,
  designObj = NULL,
  complexModel = formula(y~(Year+season)^2),
  simpleModel = formula(y~Year+season),
  indexModel = formula(y~Year+season),
  modelTry = c("Delta-Lognormal","Tweedie","Lognormal"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.05,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 1000,
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL,
  reportType = "both"
)


### original functions
# setupObj<-bycatchSetup(
#   modelTry = c("Delta-Lognormal","Tweedie"),
#   obsdat = obsdatExample,
#   logdat = logdatExample,
#   yearVar = "Year",
#   obsEffort = "sampled.sets",
#   logEffort = "sets",
#   logUnsampledEffort = NULL,
#   includeObsCatch  = FALSE,
#   matchColumn = NA,
#   factorNames = c("Year","season"),
#   randomEffects = NULL,
#   randomEffects2 = NULL,
#   EstimateIndex = FALSE,
#   EstimateBycatch = TRUE,
#   logNum = NA,
#   sampleUnit = "trips",
#   complexModel = formula(y~(Year+season)^2),
#   simpleModel = formula(y~Year),
#   indexModel = formula(y~Year),
#   designMethods = c("Ratio", "Delta"),
#   designVars = c("Year","season"),
#   designPooling = TRUE,
#   poolTypes=c("adjacent","all"),
#   pooledVar=c(NA,NA),
#   adjacentNum=c(1, NA),
#   minStrataUnit = 1,
#   baseDir = working.dir,
#   runName = "SimulatedExample",
#   runDescription = "Example with simulated data",
#   common = "Simulated species",
#   sp = "Genus species",
#   obsCatch = "Catch",
#   catchUnit = "number",
#   catchType = "discard"
# )
#
# dataCheck(setupObj)
#
# modelfit<- bycatchFit(
#   setupObj = setupObj,
#   selectCriteria = "BIC",
#   DoCrossValidation = TRUE,
#   DredgeCrossValidation = FALSE,
#   ResidualTest = FALSE,
#   CIval = 0.05,
#   VarCalc = "Simulate",
#   useParallel = TRUE,
#   nSims = 1000,
#   baseDir = getwd(),
#   plotValidation = FALSE,
#   trueVals = NULL,
#   trueCols = NULL,
#   doReport = TRUE
# )



LLSIM_BUM_Example_logbook <- LLSIM_BUM_Example_logbook[,-1]
LLSIM_BUM_Example_observer <- LLSIM_BUM_Example_observer[,-1]

shortdf_logbook <- droplevels(LLSIM_BUM_Example_logbook[LLSIM_BUM_Example_logbook$Year>2010 & LLSIM_BUM_Example_logbook$fleet==2,])
shortdf_observer <-droplevels(LLSIM_BUM_Example_observer[LLSIM_BUM_Example_observer$Year>2010 &LLSIM_BUM_Example_observer$fleet==2,])

setupObj<-bycatchSetup_new(
  obsdat = shortdf_observer,
  logdat = shortdf_logbook,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM"),
  catchUnit = "number",
  catchType = "catch",
  logNum = NA,
  sampleUnit = "trips",
  factorVariables = c("Year","area","season"),
  numericVariables = "lat",
  logUnsampledEffort = NULL,
  includeObsCatch  = FALSE,
  matchColumn = NA,
  EstimateIndex = FALSE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_3newfunc",
  runDescription = "LLSIMBUMtripExample_3newfunc",
  common = c("Swordfish","Blue marlin"),
  sp = c("Xiphias gladius","Makaira nigricans"),
  report = "html"
)


design_est <- bycatchDesign_new(
  setupObj = setupObj,
  designMethods = c("Ratio", "Delta"),
  designVars = c("Year"),
  designPooling = FALSE,
  poolTypes=c("adjacent","all"),
  pooledVar=c(NA,NA),
  adjacentNum=c(1,NA),
  minStrataUnit = 1
)


model_est <- bycatchFit_new(
  setupObj = setupObj,
  designObj = design_est,
  complexModel = formula(y~Year+season+area),
  simpleModel = formula(y~Year+season),
  indexModel = formula(y~Year+season),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.05,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 1000,
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL
)

##using original functions
setupObj<-bycatchSetup(
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
  obsdat = shortdf_observer,
  logdat = shortdf_logbook,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  factorNames = c("Year","season","area"),
  EstimateIndex = FALSE,
  EstimateBycatch = TRUE,
  logNum = NA,
  sampleUnit = "trips",
  complexModel = formula(y~Year+season+area), #upper range of models to compare with information criteria,
  simpleModel = formula(y~Year+season), #lower range of models to compare
  indexModel = formula(y~Year+season),
  designMethods =c("Ratio", "Delta"),
  baseDir = working.dir,
  runName = "LLSIMBUM_originalfunc",
  runDescription = "LLSIMBUM_originalfunc",
  common = c("Swordfish","Blue marlin")[2],
  sp = c("Xiphias gladius","Makaira nigricans")[2],
  obsCatch = c("SWO","BUM")[2],
  catchUnit = "number",
  catchType = "catch"
)

bycatchFit(
  setupObj,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  DredgeCrossValidation = FALSE,
  ResidualTest = FALSE,  #Generally should be FALSE.
  CIval = 0.05,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 100,
  baseDir = working.dir,
)




