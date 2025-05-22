# Setup required before new plot for standalone testing, not needed in markdown
# have obsdat and logdat from LLSIM
library(GGally)
factorVariables<-c("Year","season","fleet","gear")
numericVariables<-c("lat","lon","hbf")
figurenum<-1
obsdat<-LLSIM_BUM_Example_observer
logdat<-LLSIM_BUM_Example_logbook

## Code to be added to the dataCheck markdown effort comparison section
nmax<-100000  #If there are more data points than this, subsample to avoid hanging R
#Sampled 100000 if there were more than that so R won't hang
allData<-bind_rows(list(Obs=obsdat[,unique(c(factorVariables,numericVariables))],
                        Eff=logdat[,unique(c(factorVariables,numericVariables))]),
                   .id="Source")%>%arrange(Source)
if(nrow(allData)>nmax) allData<-allData[sort(sample(1:nrow(allData),nmax,replace=FALSE)),]
#Convert factors with>15 levels to number so it will plot visibly
temp<-allData %>%summarise(across(all_of(factorVariables), ~ n_distinct(.)))
temp<-names(temp)[temp[1,]>15]
allData<-allData %>%
  mutate(across(all_of(temp),~as.numeric(factor(.))))
#Make sure other factors are factors
temp<-factorVariables[!factorVariables %in% temp]
allData<-allData %>%
  mutate(across(all_of(temp),~factor(.)))
#Numerical pairs plot
cat("\n")
cat("\n Figure ",figurenum,". Effort across all pairs of numeric variables. \n",sep="")
print(ggpairs(filter(allData),
          aes(color=Source,fill=Source,alpha=0.2),
          columns=numericVariables)+
    scale_color_manual(values=c("grey","darkblue"))+
    scale_fill_manual(values=c("grey","darkblue")))
figurenum<-figurenum+1
#factor variables
cat("\n")
cat("\n Figure ",figurenum,". Effort across pairs of categorical variables. \n",sep="")
print(ggpairs(filter(allData),
        aes(color=Source,fill=Source,alpha=0.2),
        columns=factorVariables)+
  scale_color_manual(values=c("grey","darkblue"))+
  scale_fill_manual(values=c("grey","darkblue")))
figurenum<-figurenum+1
