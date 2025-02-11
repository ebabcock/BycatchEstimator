# This takes data on lengths of individual fish and converts it to 
# catches of all fish, large fish or small fish, by species
library(tidyverse)
obsdatlength<-read.csv("obsdatlength.csv")  #Made up data example
obsdatlength<-mutate(obsdatlength,Large=ifelse(length>100,1,0))  #Add a colum for 
#whether the fish is large (>100cm)
view(obsdatlength)

#The following produces set by set data, with columns
# for each species, for large fish, and for large swordfish, 
# as well as a sum of effort
obsdatInput<-obsdatlength %>% group_by(TripID,Set) %>%
  summarize(Year=min(Year),
            hooks=mean(Hooks),
            BUM=sum(Species=="BUM"),
            SWO=sum(Species=="SWO"),
            SMA=sum(Species=="SMA"),
            largeFish=sum(Large),
            largeSWO=sum(Large*(Species=="SWO")))
            
view(obsdatInput)
