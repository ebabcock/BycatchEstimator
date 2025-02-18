library(tidyverse)
library(gridExtra)
library(BycatchEstimator)
#setwd("~/Box/bycatch project (ebabcock@miami.edu)/LLSim")

#Complete dataset
allsets<-LLSIM_BUM_Example_logbook
dim(allsets)
summary(allsets)


# Make subsets for observer coverage
#Observer coverage, allocated randomly by trip
obsCoverage<-0.05

trip<-sort(unique(allsets$trip))
n<-length(trip)
tripSampled<-sample(trip,size=n*obsCoverage,replace=FALSE)

samplesets<-filter(allsets,trip %in% tripSampled )

## FOr using in tool
obsdat<-samplesets
logdat<-allsets
dim(obsdat)
dim(logdat)
