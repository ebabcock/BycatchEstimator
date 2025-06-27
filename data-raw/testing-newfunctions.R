devtools::load_all()

obsdatExample <- obsdatExample
logdatExample <- logdatExample

working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

setupObj_new <- bycatchSetup(
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
  EstimateIndex = FALSE,
  baseDir = working.dir,
  runName = "SimulatedExample_newfunctions",
  runDescription = "SimulatedExample_newfunctions",
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
  complexModel = formula(y~(Year+season)^2),
  simpleModel = formula(y~Year+season),
  indexModel = formula(y~Year+season),
  modelTry = c("Delta-Lognormal","Tweedie","Lognormal"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.05,
  VarCalc = "DeltaMethod",
  useParallel = TRUE,
  nSims = 100,
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL,
  reportType = "html"
)

### LLSIM data

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
  nSims = 100,
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL
)
