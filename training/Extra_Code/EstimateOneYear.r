library(BycatchEstimator)
library(MuMIn)

#This code shows how to use the model with only one year of 
# data.  You still need a Year variable becuase some 
# elements of the code expect it. Year does not have to 
# be one of the model variables, but it must be in 
# designVars if you are using the desgin based estimators

obsdatOne<-filter(LLSIM_BUM_Example_observer,Year==2000 )
logdatOne<-filter(LLSIM_BUM_Example_logbook,Year==2000 )
table(obsdatOne$month,obsdatOne$area)

setupObj<-bycatchSetup(
  modelTry = c("TMBdelta-Lognormal"),
  obsdat = obsdatOne,
  logdat = logdatOne,
  yearVar = "Year",
  obsEffort = "hooks",
  logEffort = "hooks",
  factorNames = c("month","area","fleet"), #To be sure that factors are not interpreted as numbers
  EstimateIndex = FALSE,
  EstimateBycatch = TRUE,  #Should be true
  logNum = NA,
  sampleUnit = "trips",
  complexModel = formula(y~month+area+fleet),  #Complex model not used in this set
  simpleModel = formula(y~month), # This defines strata for the models, and also some output files
  indexModel = NA,  # not used
  designMethods =c("Ratio","Delta"),  #THis specifies the design-based methods to use
  designVars=c("Year","month","area"),  #Stratification variables
  designPooling = FALSE,  #No pooling.  #Bycatch will be estimated at 0 for strata with no observer data
  baseDir = getwd(),
  runName = "LLSIMBUMOneYear",
  runDescription = "LLSIm BUM for 2000 only",
  common = "Blue marlin",
  sp = "Makaira nigricans",
  obsCatch = "BUM",
  catchUnit = "number",
  catchType = "catch"
)
dataCheck(setupObj)
bycatchFit(setupObj)

#If we wanted results by month (or any other variable) within the year
# these are not printed automatically, but they are in the .csv outputs
# if month was included in simpleModel. To plot

nbMonth<-read.csv("Output LLSIMBUMOneYear/Blue marlin catch/Blue marlincatchTMBdelta-LognormalStratumSummary.csv")
nbMonth<-mutate(nbMonth,Method="Neg Binomial")%>%
  select(month=strata,Mean=Total,SE=Total.se,LCI=TotalLCI,UCI=TotalUCI) %>%
  mutate(Method="NegBin")
view(nbMonth)

designStrata<-read.csv("Output LLSIMBUMOneYear/Blue marlin catch/Blue marlincatchDesignStrata.csv") 
view(designStrata)
#THis is by strata including month and area, so we need to do some adding

designMonth<-designStrata %>% group_by(month)%>%
  summarize(ratioMean=sum(ratioMean),
           deltaMean=sum(deltaMean),
           ratioSE=sqrt(sum(ratioSE^2)),
           deltaSE=sqrt(sum(deltaSE^2)))
view(designMonth)

yearSumLong<-pivot_longer(designMonth,ratioMean:deltaSE,names_to = "Method") %>%
  mutate(parameter=ifelse(grepl("Mean",Method),"Mean","SE"),
         Method=sub("Mean","",Method),
         Method=sub("SE","",Method)) %>%
  pivot_wider(names_from=parameter,values_from=value) %>%
  mutate(LCI=Mean+2*SE,
         UCI=Mean-2*SE)%>%
  mutate(UCI=ifelse(UCI>0,UCI,0)) %>%
  bind_rows(nbMonth)

ggplot(yearSumLong,
       aes(x=month,y=Mean,ymin=Mean-SE,ymax=Mean+SE,color=Method,fill=Method))+
  geom_line()+
  geom_point()+
  geom_ribbon(alpha=0.5)
