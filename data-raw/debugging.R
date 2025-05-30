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
  VarCalc = "None",
  useParallel = TRUE,
  nSims = 1000,
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL,
  reportType = "html"
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


### matching trips
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest,
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
  logUnsampledEffort = "unsampledEffort",
  includeObsCatch  = TRUE,
  matchColumn = "trip",
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_tripmatching&index",
  runDescription = "LLSIMBUMtripExample_tripmatching&index",
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
  complexModel = formula(y~Year+area),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year+area),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.05,
  VarCalc = "DeltaMethod",
  useParallel = TRUE,
  nSims = 1000,
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL
)


### validation data
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)

trueVals<-read.csv(paste0(working.dir,"TotalAnnualCatches.csv"))
trueVals<-filter(trueVals,Year %in% 2011:2015)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[1],
  catchUnit = "number",
  catchType = "catch",
  logNum = NA,
  sampleUnit = "trips",
  factorVariables = c("Year","area","fleet"),
  numericVariables = "lat",
  logUnsampledEffort = "unsampledEffort",
  includeObsCatch  = TRUE,
  matchColumn = "trip",
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_validation&novar",
  runDescription = "LLSIMBUMtripExample_validation&novar",
  common = c("Swordfish","Blue marlin")[1],
  sp = c("Xiphias gladius","Makaira nigricans")[1],
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
  complexModel = formula(y~Year+area+fleet),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year+area),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie","TMBnbinom2"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.05,
  VarCalc = "None",
  useParallel = TRUE,
  nSims = 1000,
  plotValidation = TRUE,
  trueVals = trueVals,
  trueCols = c("Total.BUM","Total.SWO")[2]
)


### cpue random effects
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2001:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2001:2015)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[2],
  catchUnit = "number",
  catchType = "catch",
  logNum = NA,
  sampleUnit = "trips",
  factorVariables = c("Year","area","fleet","season"),
  numericVariables = NA,
  logUnsampledEffort = NULL,
  includeObsCatch  = FALSE,
  matchColumn = NA,
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_cpue_randomeff",
  runDescription = "LLSIMBUMtripExample_cpue_randomeff",
  common = c("Swordfish","Blue marlin")[2],
  sp = c("Xiphias gladius","Makaira nigricans")[2],
  report = "html"
)

model_est <- bycatchFit_new(
  setupObj = setupObj,
  complexModel = formula(y~Year+area+fleet+season+area:season),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year+area),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie","TMBnbinom2","TMBdelta-Gamma"),
  randomEffects = c("Year:area"),
  randomEffects2 = c("Year:area"),
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.1,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 100,
  plotValidation = FALSE
)

### cpue add more variables
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2001:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2001:2015)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[2],
  catchUnit = "number",
  catchType = "catch",
  logNum = NA,
  sampleUnit = "trips",
  factorVariables = c("Year","area","fleet","season","light"),
  numericVariables = c("hbf","habBUM"),
  logUnsampledEffort = NULL,
  includeObsCatch  = FALSE,
  matchColumn = NA,
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_cpue_morevariab",
  runDescription = "LLSIMBUMtripExample_cpue_morevariab",
  common = c("Swordfish","Blue marlin")[2],
  sp = c("Xiphias gladius","Makaira nigricans")[2],
  report = "html"
)

model_est <- bycatchFit_new(
  setupObj = setupObj,
  complexModel = formula(y~Year+area+fleet+season+light+hbf+habBUM),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year+area),
  modelTry = c("TMBdelta-Lognormal","TMBnbinom1","TMBdelta-Gamma","TMBnbinom2","TMBtweedie"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.05,
  VarCalc = "DeltaMethod",
  useParallel = TRUE,
  nSims = 100,
  plotValidation = FALSE
)

### cpue no bycatch
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2001:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2001:2015)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = NULL,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[2],
  catchUnit = "number",
  catchType = "catch",
  logNum = NA,
  sampleUnit = "trips",
  factorVariables = c("Year","area","fleet","season"),
  numericVariables = NA,
  logUnsampledEffort = NULL,
  includeObsCatch  = FALSE,
  matchColumn = NA,
  EstimateIndex = TRUE,
  EstimateBycatch = FALSE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_cpue_nobycatch",
  runDescription = "LLSIMBUMtripExample_cpue_nobycatch",
  common = c("Swordfish","Blue marlin")[2],
  sp = c("Xiphias gladius","Makaira nigricans")[2],
  report = "html"
)

model_est <- bycatchFit_new(
  setupObj = setupObj,
  complexModel = formula(y~Year+area+fleet+season),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year+area),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie","TMBnbinom2","TMBdelta-Gamma"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.1,
  VarCalc = "Simulate",
  useParallel = TRUE,
  nSims = 100,
  plotValidation = FALSE
)


### matching trips - observer trips not matched
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)
#adding some fake trips to observer data
fake.trips <- obsdatTest[5:6,]
fake.trips$trip[1] <- "22229999988888"
fake.trips$trip[2] <- "22229999955555"
obsdatTest2 <- rbind(fake.trips,obsdatTest)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest2,
  logdat = logdatTest,
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
  logUnsampledEffort = "unsampledEffort",
  includeObsCatch  = TRUE,
  matchColumn = "trip",
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_tripmatchingerror",
  runDescription = "LLSIMBUMtripExample_tripmatchingerror",
  common = c("Swordfish","Blue marlin"),
  sp = c("Xiphias gladius","Makaira nigricans"),
  report = "html"
)

# design_est <- bycatchDesign_new(
#   setupObj = setupObj,
#   designMethods = c("Ratio", "Delta"),
#   designVars = c("Year"),
#   designPooling = FALSE,
#   poolTypes=c("adjacent","all"),
#   pooledVar=c(NA,NA),
#   adjacentNum=c(1,NA),
#   minStrataUnit = 1
# )
#
# model_est <- bycatchFit_new(
#   setupObj = setupObj,
#   complexModel = formula(y~Year+area),
#   simpleModel = formula(y~Year),
#   indexModel = formula(y~Year+area),
#   modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
#   randomEffects = NULL,
#   randomEffects2 = NULL,
#   selectCriteria = "BIC",
#   DoCrossValidation = TRUE,
#   CIval = 0.05,
#   VarCalc = "DeltaMethod",
#   useParallel = TRUE,
#   nSims = 1000,
#   plotValidation = FALSE,
#   trueVals = NULL,
#   trueCols = NULL
# )


### matching trips - observer trips with NAs
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)
#changing some trips to NAs
natrips <- obsdatTest[5:10,]
natrips$SWO <- NA
obsdatTest2 <- rbind(obsdatTest,natrips)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest2,
  logdat = logdatTest,
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
  logUnsampledEffort = "unsampledEffort",
  includeObsCatch  = TRUE,
  matchColumn = "trip",
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_tripmatchingNA",
  runDescription = "LLSIMBUMtripExample_tripmatchingNA",
  common = c("Swordfish","Blue marlin"),
  sp = c("Xiphias gladius","Makaira nigricans"),
  report = "html"
)

### matching trips - logbook trips with NAs
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)
#changing some trips to NAs
natrips <- logdatTest[5:10,]
natrips$SWO <-NA
logdatTest2 <- rbind(logdatTest,natrips)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest2,
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
  logUnsampledEffort = "unsampledEffort",
  includeObsCatch  = TRUE,
  matchColumn = "trip",
  EstimateIndex = TRUE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_logbookNA",
  runDescription = "LLSIMBUMtripExample_logbookNA",
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
  complexModel = formula(y~Year+area),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year+area),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
  randomEffects = NULL,
  randomEffects2 = NULL,
  selectCriteria = "BIC",
  DoCrossValidation = TRUE,
  CIval = 0.05,
  VarCalc = "DeltaMethod",
  useParallel = TRUE,
  nSims = 1000,
  plotValidation = FALSE,
  trueVals = NULL,
  trueCols = NULL,
  reportType = "html"
)


### missing levels in observer data
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)
#exclude one season and area from observer data
obsdatTest2 <- filter(obsdatTest, season %in% c(1,2,3) & area == "N")

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest2,
  logdat = logdatTest,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[1],
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
  runName = "LLSIMBUMtripExample_missingobserverdata",
  runDescription = "LLSIMBUMtripExample_missingobserverdata",
  common = c("Swordfish","Blue marlin")[1],
  sp = c("Xiphias gladius","Makaira nigricans")[1],
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
  complexModel = formula(y~Year+area+season),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
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

### missing levels in logbook data
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)
#exclude one season and area from logbook data
logdatTest2 <- filter(logdatTest, season %in% c(1,2,3) & area == "N")

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest2,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[1],
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
  runName = "LLSIMBUMtripExample_missinglogdata",
  runDescription = "LLSIMBUMtripExample_missinglogdata",
  common = c("Swordfish","Blue marlin")[1],
  sp = c("Xiphias gladius","Makaira nigricans")[1],
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
  complexModel = formula(y~Year+area+season),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
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


### NAs in observer data, not matching trips
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)
#one season and area from observer data as NA
obsdatTest2 <- obsdatTest
obsdatTest2$season[obsdatTest2$season==4] <- NA
obsdatTest2$area[obsdatTest2$area=="S"] <- NA

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest2,
  logdat = logdatTest,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[1],
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
  runName = "LLSIMBUMtripExample_missingobserverdataNA",
  runDescription = "LLSIMBUMtripExample_missingobserverdataNA",
  common = c("Swordfish","Blue marlin")[1],
  sp = c("Xiphias gladius","Makaira nigricans")[1],
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
  complexModel = formula(y~Year+area+season),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
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


## testing adding R-squared in findBestModelFunc
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[1],
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
  runName = "LLSIMBUMtripExample_Rsquared",
  runDescription = "LLSIMBUMtripExample_Rsquared",
  common = c("Swordfish","Blue marlin")[1],
  sp = c("Xiphias gladius","Makaira nigricans")[1],
  report = "html"
)

model_est <- bycatchFit_new(
  setupObj = setupObj,
  complexModel = formula(y~Year+area+season),
  simpleModel = formula(y~Year),
  indexModel = formula(y~Year),
  modelTry = c("TMBdelta-Lognormal","Delta-Lognormal","TMBtweedie"),
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


### testing NA warning messages
working.dir <- "C:/Users/aadao/OneDrive/Documents/NA work 2025/ICCAT work/Tool improvements project/"

obsdatTest<-filter(LLSIM_BUM_Example_observer,Year %in% 2011:2015)
logdatTest<-filter(LLSIM_BUM_Example_logbook,Year %in% 2011:2015)
#changing some trips to NAs
natrips <- obsdatTest[5:10,]
natrips$SWO <- natrips$area<- NA
obsdatTest2 <- rbind(obsdatTest,natrips)
natrips<-logdatTest[5:10,]
natrips$SWO <- natrips$area<- NA
logdatTest2 <- rbind(logdatTest,natrips)

setupObj<-bycatchSetup_new(
  obsdat = obsdatTest,
  logdat = logdatTest2,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM")[1],
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
  runName = "LLSIMBUMtripExample_NAmessage",
  runDescription = "LLSIMBUMtripExample_NAmessage",
  common = c("Swordfish","Blue marlin")[1],
  sp = c("Xiphias gladius","Makaira nigricans")[1],
  report = "html"
)
