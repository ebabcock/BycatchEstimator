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
  runName = "SimulatedExample_new2functions",
  runDescription = "SimulatedExample_new2functions",
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



LLSIM_BUM_Example_logbook <- LLSIM_BUM_Example_logbook[,-1]
LLSIM_BUM_Example_observer <- LLSIM_BUM_Example_observer[,-1]

setupObj<-bycatchSetup_new(
  obsdat = droplevels(LLSIM_BUM_Example_observer[LLSIM_BUM_Example_observer$Year>2010 &LLSIM_BUM_Example_observer$fleet==2,]),
  logdat = droplevels(LLSIM_BUM_Example_logbook[LLSIM_BUM_Example_logbook$Year>2010 & LLSIM_BUM_Example_logbook$fleet==2,]),
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  obsCatch = c("SWO","BUM"),
  catchUnit = "number",
  catchType = "catch",
  logNum = NA,
  sampleUnit = "trips",
  factorVariables = c("Year","area"),
  numericVariables = "lat",
  logUnsampledEffort = NULL,
  includeObsCatch  = FALSE,
  matchColumn = NA,
  EstimateIndex = FALSE,
  EstimateBycatch = TRUE,
  baseDir = working.dir,
  runName = "LLSIMBUMtripExample_2newfunc",
  runDescription = "LLSIMBUMtripExample_2newfunc",
  common = c("Swordfish","Blue marlin"),
  sp = c("Xiphias gladius","Makaira nigricans"),
  report = "html"
)


design_est <- bycatchDesign_new(
  setupObj = setupObj,
  designMethods = c("Ratio"),
  designVars = c("Year"),
  designPooling = FALSE,
  poolTypes=c("adjacent","all"),
  pooledVar=c(NA,NA),
  adjacentNum=c(1,NA),
  minStrataUnit = 1
)


setupDesign_new<-readRDS(file=paste0(working.dir,
                                     "Output LLSIMBUMtripExample_2newfunc/","2025-05-06_BycatchDesignSpecification.rds"))

