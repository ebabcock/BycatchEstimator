
#' Basic ratio estimator with variance (Cochran)
#'
#' Output is mean and standard error of bycatch by stratum and the total bycatch with SE. Assumes unobserved strata have zero catch
#' @param x x, y and g are vectors giving the effort/catch, bycatch and stratum of each observed sample unit.
#' @param y x, y and g are vectors giving the effort/catch, bycatch and stratum of each observed sample unit.
#' @param g x, y and g are vectors giving the effort/catch, bycatch and stratum of each observed sample unit.
#' @param X X is the total effort/catch by stratum
#' @param N N is the total number of sample units by stratum, if available, otherwise total effort
#' @param G Value
#' @importFrom stats var cov
#' @keywords internal
ratio.func= function(x,y,g,X,N,G) {
  N=N[G %in% g]
  X=X[G %in% g]
  G=G[G %in% g]
  a=order(G)
  X=X[a]
  N=N[a]
  G=G[a]
  if(sum(X-N)!=0)  n=tapply(x,g,length)  else  n=tapply(x,g,sum)
  f=n/N
  xmean=tapply(x,g,mean)
  ymean=tapply(y,g,mean)
  Rhat=ymean/xmean
  sx2=var(x)
  sy2=var(y)
  sxy=cov(x,y)
  stratum.est=X*Rhat
  stratum.var=X^2*(1-f)/(n*xmean^2)*(sy2+Rhat^2*sx2-2*Rhat*sxy)
  total.est=sum(stratum.est)
  total.var=sum(stratum.var)
  list(stratum.est=stratum.est,stratum.se=sqrt(stratum.var),coverage=f,total.est=total.est,total.se=sqrt(total.var))
}

#' Function for Pennington(1983) method
#'
#' Output is result of Gm(t) function for use in calculations of the design based delta estimator
#' @param m m is the number of non-zero values
#' @param t t is a real number input
#' @param jmax jmax is the upper limit of an infinute sum used in calculations. 10 is usually sufficient
#' @keywords internal
Gm<-function(m,t,jmax=10) {
  x=c(NA,m+2*(2:jmax)-3)
  func1=function(j) prod(seq(m+1,x[j],by=2))
  y=c(NA,sapply(2:jmax,func1))
  j=2:jmax
  1+(m-1)*t/m+sum((m-1)^(2*j-1)*(t/m)^j/(factorial(j)*y[j]))
}

#' Function for Pennington(1983) method mean
#'
#' Output is mean of the delta estimator
#' @param x is vector of data input
#' @importFrom stats var
#' @keywords internal
deltaEstimatorMean<-function(x) {
  x=x[!is.na(x)]
  n=length(x)
  r=length(x[x==0])
  p=(n-r)/n
  ybar=mean(log(x[x>0]))
  s2=var(log(x[x>0]))
  returnval=p*exp(ybar)*Gm(m=n-r,t=s2/2)
  if(r==n-1) returnval=x[x>0]/n
  if(r==n) returnval=0
  if(n==0) returnval=NA
  returnval
}

#' Function for Pennington(1983) method variance
#'
#' Output is variance of the delta estimator
#' @param x is vector of data input
#' @importFrom stats var
#' @keywords internal
deltaEstimatorVar<-function(x) {
  x=x[!is.na(x)]
  n=length(x)
  r=length(x[x==0])
  p=(n-r)/n
  ybar=mean(log(x[x>0]))
  s2=var(log(x[x>0]))
  returnval=p*exp(2*ybar)*(Gm(m=n-r,t=2*s2)-(n-r-1)/(n-1)*
                             Gm(m=n-r,t=s2*(n-r-2)/(n-r-1)))
  if(r==n-1) returnval=(x[x>0])^2/n
  if(r==n) returnval=0
  if(n<2) returnval=NA
  if(n==0) returnval=NA
  returnval
}

#' Function for Pennington(1983) method SE of the mean squAred
#'
#' Output is SE squared of the delta estimator
#' @param x is vector of data input
#' @importFrom stats var
#' @keywords internal
deltaEstimatorSE2<-function(x) {
  x=x[!is.na(x)]
  n=length(x)
  r=length(x[x==0])
  p=(n-r)/n
  ybar=mean(log(x[x>0]))
  s2=var(log(x[x>0]))
  returnval=p*exp(2*ybar)*(p*Gm(m=n-r,t=s2/2)^2-(n-r-1)/(n-1)*
                             Gm(m=n-r,t=s2*(n-r-2)/(n-r-1)))
  if(r==n-1) returnval=(x[x>0]/n)^2
  if(r==n) returnval=0
  if(n<2) returnval=NA
  if(n==0) returnval=NA
  returnval
}




#'Function to find mode of a categorical variable
#'
#' @param x value
#' @importFrom stats aggregate
#' @keywords internal
mostfreqfunc<-function(x) {
  x=x[!is.na(x)]
  if(length(x)>0) {
    y=aggregate(x,list(x),length)
    temp=y$Group.1[y$x==max(y$x)][1]
  } else temp=NA
  temp
}

#'Function to count outliers, defined as more than numSD standard deviations from the mean.
#'
#' @param x Value
#' @param numSD Value
#' @importFrom stats sd
#' @keywords internal
outlierCountFunc=function(x,numSD=8) {
  length(which(x>mean(x,na.rm=TRUE)+numSD*sd(x,na.rm=TRUE) |  x<mean(x,na.rm=TRUE)-numSD*sd(x,na.rm=TRUE)))
}

#'Standard error of a mean
#'
#' @param x Value
#' @importFrom stats sd
#' @keywords internal
standard.error<-function(x) {
  x=x[!is.na(x)]
  sd(x)/sqrt(length(x))
}

#' Function to find best model by information criteria, by model type
#'
#' @param obsdatval Value
#' @param modType Value
#' @param requiredVarNames Value
#' @param allVarNames Value
#' @param complexModel Value
#' @param randomEffects Value
#' @param useParallel Value
#' @param selectCriteria Value
#' @param varExclude Value
#' @param printOutput Value
#' @param catchType Value
#' @param common Value
#' @param dirname Value
#' @param run Value
#' @param modelScenario Value
#' @import MuMIn parallel doParallel tweedie glmmTMB
#' @importFrom reshape2 colsplit
#' @importFrom stats anova na.fail as.formula coef Gamma glm.control formula lm glm vcov
#' @importFrom MASS glm.nb
#' @keywords internal
findBestModelFunc<-function(obsdatval, modType, requiredVarNames, allVarNames, complexModel,
  randomEffects=NULL, useParallel, selectCriteria, varExclude, printOutput=FALSE,
  catchType = NULL, common = NULL, dirname = NULL, run = NULL,modelScenario=NULL) {

  offset<-TMBfamily<-NULL
  requiredVarNames<-requiredVarNames[!requiredVarNames %in% varExclude]
  if(length(requiredVarNames)>0) keepVars=requiredVarNames else keepVars=NULL
  extras=c("AICc","AIC", "BIC")
  if(!is.null(randomEffects)) randomEffects<-paste0("(1|",randomEffects,")")
  #Check if parallel is possible
  if(useParallel){
    NumCores<-detectCores()  #Check if machine has multiple cores for parallel processing
    if(NumCores >= 3){
      cl2<-makeCluster(NumCores-2)
      registerDoParallel(cl2)
      #For MuMIn
      #cl2<-makeCluster(NumCores-2)
      # clusterEvalQ(cl2, {
      #   library(glmmTMB)
      #   library(cplm)
      #   library(MASS)
      # })
    } else {
      useParallel <- FALSE
    }
  }
  #Set up input data
  if(modType %in% c("Binomial","TMBbinomial"))
    obsdatval$y=ifelse(obsdatval$Catch>0,1,0)
  if(modType %in% c("NegBin","TMBnbinom1","TMBnbinom2"))
    obsdatval$y=round(obsdatval$Catch)
  if(modType %in% c("Tweedie","TMBtweedie","Normal","TMBnormal"))
    obsdatval$y=obsdatval$cpue
  if(modType %in% c("Lognormal","TMBlognormal"))
    obsdatval$y=log(obsdatval$cpue+0.1)
  if(modType %in% c("Gamma","TMBgamma") )
    obsdatval$y=obsdatval$cpue+0.1
  if(modType %in% c("Delta-Lognormal","TMBdelta-Lognormal")) {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
  }
  if(modType %in% c("Delta-Gamma","TMBdelta-Gamma")) {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
  }
  funcName=case_when(modType %in% c("Normal","Lognormal","Delta-Lognormal")~"lm",
                     modType %in% c("Binomial","Gamma","Delta-Gamma")~"glm",
                     modType %in% c("NegBin")~"glm.nb",
                     modType %in% c("Tweedie") ~"cpglm",
                     grepl("TMB",modType)~"glmmTMB")
  args=list(formula="",data=obsdatval,na.action=na.fail)
  if(funcName=="glm") args=c(args,list(control=list(epsilon = 1e-6,maxit=100)))
  if(modType %in% c("Binomial","TMBbinomial")) {
    args=c(args,list(family="binomial"))
  }
  if(modType %in% c("Gamma","TMBgamma","Delta-Gamma","TMBdelta-Gamma"))  {
    args=c(args,list(family=Gamma(link="log")))
    TMBfamily=Gamma(link="log")
  }
  if(modType %in% c("NegBin")) {
    offset="+offset(log(Effort))"
    keepVars=c(requiredVarNames,"offset(log(Effort))")
  }
  if(modType %in% c("TMBnbinom1","TMBnbinom2") ){
    TMBfamily=gsub("TMB","",modType)
    offset="+offset(log(Effort))"
    keepVars=paste0("cond(",c(requiredVarNames,"offset(log(Effort))"),")")
    args=c(args,family=TMBfamily)
  }
  if(modType %in% c("TMBtweedie") ){
    TMBfamily=gsub("TMB","",modType)
    if(length(requiredVarNames)>0) keepVars=paste0("cond(",requiredVarNames,")")
    args=c(args,list(family=TMBfamily))
  }
  if(modType %in% c("TMBnormal","TMBlognormal","TMBdelta-Lognormal")) {
    TMBfamily="gaussian"
    if(length(requiredVarNames)>0) keepVars=paste0("cond(",requiredVarNames,")")
  }
  if(modType %in% c("TMBbinomial")){
    TMBfamily="binomial"
    if(length(requiredVarNames)>0) keepVars=paste0("cond(",requiredVarNames,")")
  }
  if(modType %in% c("TMBgamma","TMBdelta-Gamma")) {
    if(length(requiredVarNames)>0) keepVars=paste0("cond(",requiredVarNames,")")
  }
  allVarNames<-as.vector(getAllTerms(complexModel))
  allVarNames<-allVarNames[!allVarNames %in% varExclude]
  if(length(requiredVarNames)>0 )
    formulaList<-list(as.formula(paste("y~",paste(c(allVarNames,randomEffects),collapse="+"),offset)),
                    as.formula(paste("y~",paste(c(allVarNames[!grepl(":",allVarNames)],randomEffects),collapse="+"),offset)),
                    as.formula(paste("y~",paste(c(allVarNames[!grepl(":",allVarNames) &!allVarNames %in% varExclude],randomEffects),collapse="+"),offset)),
                    as.formula(paste("y~",paste(requiredVarNames,collapse="+"),offset)),NA) else
        formulaList<-list(as.formula(paste("y~",paste(c(allVarNames,randomEffects),collapse="+"),offset)),
                    as.formula(paste("y~",paste(c(allVarNames[!grepl(":",allVarNames)],randomEffects),collapse="+"),offset)),
                    as.formula(paste("y~",paste(c(allVarNames[!grepl(":",allVarNames) &!allVarNames %in% varExclude],randomEffects),collapse="+"),offset)),
                    NA,NA)
  args$formula=formulaList[[1]]
  if(! modType=="Tweedie") modfit1<-try(do.call(funcName,args))  else
    modfit1<-try(cplm::cpglm(formulaList[[1]],data=obsdatval,na.action=na.fail))
  for(i in 2:(length(formulaList))-1) {
    if(class(modfit1)[1] %in% c("glm","lm","glm.nb")) {
      if(modfit1$rank<length(coef(modfit1))) class(modfit1)<-"try-error"
    }
    if(class(modfit1)[1] == "cpglm")  {
      if(length(coef(modfit1))!=dim(modfit1$vcov)[1])  class(modfit1)<-"try-error"
    }
    if(class(modfit1)[1]=="try-error")   {
      args$formula=formulaList[[i]]
      if(! modType=="Tweedie")
        modfit1<-try(do.call(funcName,args)) else
          modfit1<-try(cplm::cpglm(formulaList[[i]],data=obsdatval,na.action=na.fail))
    }
  }
  if(class(modfit1)[1]=="try-error")   {
    returnval=NULL
    print(paste(common[run],modType,"failed to converge"))
  } else {
    if(modType=="Binomial") modfit1<-glm(formula(modfit1),data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=100),na.action=na.fail)
    if(modType %in% c("Normal","Lognormal","Delta-Lognormal")) modfit1<-lm(formula(modfit1),data=obsdatval,na.action=na.fail)
    if(modType %in% c("Gamma","Delta-Gamma")) modfit1<-glm(formula(modfit1),data=obsdatval,family=Gamma(link="log"),na.action=na.fail)
    if(modType=="NegBin") modfit1<-glm.nb(formula(modfit1),data=obsdatval,control=glm.control(epsilon=1E-6,maxit=30),na.action=na.fail)
    if(modType=="Tweedie") modfit1<-cplm::cpglm(formula(modfit1),data=obsdatval,na.action=na.fail)
    if(grepl("TMB",modType) )
      modfit1<-glmmTMB(formula(modfit1),family=TMBfamily,data=obsdatval,na.action=na.fail)
    if(useParallel) {
      clusterEvalQ(cl2, {rm(list=ls())} )
      clusterExport(cl2,c("obsdatval","modfit1","keepVars","extras","TMBfamily","offset","selectCriteria"),envir=environment())
      modfit2<-MuMIn:::.dredge.par(modfit1,rank=selectCriteria,fixed=keepVars,extra=extras,cluster=cl2)
      stopCluster(cl2)
    } else {
      modfit2<-try(dredge(modfit1,rank=selectCriteria,fixed=keepVars,extra=extras))
    }
    if(class(modfit2)[1]!="try-error") {
      modfit3<-get.models(modfit2,1)[[1]]
    } else {
      modfit2<-NULL
      modfit3<-NULL
    }
    selTable<-data.frame(modfit2)
    if(!is.null(modfit2)) {
      if(is.na(sum(modfit2$weight))) {
        selTable$delta<-selTable$selectCriteria-min(selTable$selectCriteria,na.rm=TRUE)
        selTable$weight<-exp(-0.5*selTable$delta)/sum(exp(-0.5*selTable$delta),na.rm=TRUE)
      }
      selTable$R2<-addR2(modfit2,obsdatval,funcName)
    }
    if(printOutput & !is.null(modfit2)) {
      write.csv(selTable,paste0(dirname[[run]],common[run],catchType[run],modelScenario,"ModelSelection",modType,".csv"), row.names = FALSE)
      #if(modType %in% c("Binomial","NegBin")) anova1=anova(modfit3,test="Chi")
      #if(modType =="Tweedie" | grepl("TMB",modType)) anova1=NULL
      #if(modType %in% c("Normal","Lognormal","Gamma","Delta-Lognormal","Delta-Gamma")) anova1=anova(modfit3,test="F")
      # if(!is.null(anova1)) {
      #   write.csv(anova1,paste0(dirname[[run]],common[run],catchType[run],modType,"Anova.csv"), row.names = FALSE)
      # }
    }
    returnval=list(modfit3,selTable)
  }
  returnval
}


#' Generate standard errors and confidence intervals of predictions from simulation from regression coefficients and their var/covar matrix
#'
#' @param modfit1 Value
#' @param modfit2 Value
#' @param newdat Value
#' @param modtype Value
#' @param obsdatval Value
#' @param includeObsCatch Value
#' @param nsim Value
#' @param requiredVarNames Value
#' @param CIval Value
#' @param printOutput Value
#' @param catchType Value
#' @param common Value
#' @param dirname Value
#' @param run Value
#' @param randomEffects Value
#' @param randomEffects2 Value
#' @param modelScenario Value
#' @import tidyr
#' @importFrom stats predict model.matrix rbinom sigma rnorm rlnorm rnbinom quantile
#' @importFrom MASS mvrnorm gamma.shape
#' @keywords internal
makePredictionsSimVarBig<-function(modfit1, modfit2=NULL, newdat, modtype, obsdatval,
                                   includeObsCatch, nsim, requiredVarNames, CIval, printOutput=TRUE,
                                   catchType, common, dirname, run,randomEffects,randomEffects2,modelScenario) {
  #Separate out sample units
  if(includeObsCatch)    newdat$Effort=newdat$unsampledEffort/newdat$SampleUnits else
    newdat$Effort=newdat$Effort/newdat$SampleUnits
  newdat=uncount(newdat,.data$SampleUnits)
  newdatall=newdat
  #Set up output dataframes
  years=sort(unique(newdat$Year))
  yearpred=expand.grid(Year=years,Total=NA,TotalVar=NA,Total.mean=NA,TotalLCI=NA,TotalUCI=NA,Total.se=NA,Total.cv=NA)
  stratapred=expand.grid(strata=unique(newdatall$strata),Total=0,TotalVar=0,Total.mean=0,TotalLCI=NA,TotalUCI=NA,Total.se=NA,Total.cv=NA)
  stratapred$Year=newdatall$Year[match(stratapred$strata,newdatall$strata)]
  for(i in 1:length(years)) {
    newdat = newdatall[newdatall$Year==years[i],]
    nObs= nrow(newdat)
    #Get predictions
    if(modtype=="Tweedie" ) response1<-data.frame(cplm::predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(grepl("TMB",modtype)) response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE,allow.new.levels=TRUE))
    if(!grepl("TMB",modtype) & !modtype=="Tweedie") response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(dim(response1)[2]==1) {
      names(response1)="fit"
      if(modtype=="Tweedie")
        response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp",offsetval=NULL, nsim=nsim) else
          response1$se.fit=rep(NA,dim(response1)[2])
    }
    if(!is.null(modfit2))  {
      if(grepl("TMB",modtype)) response2<-data.frame(predict(modfit2,newdata=newdat,type="response",se.fit=TRUE,allow.new.levels=TRUE))
      if(!grepl("TMB",modtype) ) response2<-data.frame(predict(modfit2,newdata=newdat,type="response",se.fit=TRUE))
      names(response2)=paste0(names(response2),"2")
    }
    if(!any(is.na(response1$se.fit)) & !max(response1$se.fit/response1$fit,na.rm=TRUE)>10000)  {
      #Set up model matrices for simulation
      yvar=sub( " ", " ",formula(modfit1) )[2]
      newdat<-cbind(y=rep(1,nObs),newdat)
      names(newdat)[1]=yvar
      a=model.matrix(formula(modfit1,fixed.only=TRUE),data=newdat)
      if(!is.null(modfit2)) {
        yvar=sub( " ", " ",formula(modfit2) )[2]
        if(! yvar %in% names(newdat)) {
          newdat<-cbind(y=rep(1,nObs),newdat)
          names(newdat)[1]=yvar
        }
        b=model.matrix(formula(modfit2,fixed.only=TRUE),data=newdat)
      }
      #Get predictions sim
      if(grepl("TMB",modtype)) {
        coefvals1<-fixef(modfit1)[[1]]
        vcovvals1<-vcov(modfit1)[[1]]
      }  else
        if(!modtype=="Tweedie") {
          coefvals1<-coef(modfit1)
          vcovvals1<-vcov(modfit1)
        }
      if(!is.null(modfit2)) {
        if(grepl("TMB",modtype)) {
          coefvals2<-fixef(modfit2)[[1]]
          vcovvals2<-vcov(modfit2)[[1]]
        }  else  {
          coefvals2<-coef(modfit2)
          vcovvals2<-vcov(modfit2)
        }
      }
      NewRandomVals<-rep(0,nrow(newdat))
      NewRandomVals2<-rep(0,nrow(newdat))
      if(!is.null(randomEffects))  {
        for(j in 1:length(randomEffects)) {
          RandomVals1<-ranef(modfit1)[[1]][[j]]
          x<-match(newdat[,randomEffects[j]],rownames(RandomVals1))
          NewRandomVals[!is.na(x)]<-NewRandomVals[!is.na(x)]+RandomVals1[x[!is.na(x)],1]
        }
        if(!is.null(randomEffects2) & !is.null(modfit2))
          for(j in 1:length(randomEffects2)) {
            RandomVals2<-ranef(modfit2)[[1]][[j]]
            x<-match(newdat[,randomEffects2[j]],rownames(RandomVals2))
            NewRandomVals2[!is.na(x)]<-NewRandomVals2[!is.na(x)]+RandomVals2[x[!is.na(x)],1]
          }
      }
      if(modtype %in% c("Binomial","TMBbinomial") ){
        allpred<-cbind(newdat,response1) %>%
          mutate(Total=.data$fit,TotalVar=.data$se.fit^2+.data$fit*(1-.data$fit))
        sim=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coefvals1,vcovvals1)+NewRandomVals ) ) )
      }
      if(modtype %in% c("Normal","TMBnormal")) {
        allpred<-cbind(newdat,response1)   %>%
          mutate(Total=.data$Effort*.data$fit,
                 TotalVar=.data$Effort^2*(.data$se.fit^2+sigma(modfit1)^2))
        sim=replicate(nsim,rnorm(nObs,
                                 mean=as.vector(a %*% mvrnorm(1,coefvals1,vcovvals1)+NewRandomVals),
                                 sd=sigma(modfit1)))*newdat$Effort
      }
      if(modtype %in% c("Lognormal","TMBlognormal") ){
        allpred<-cbind(newdat,response1)   %>%
          mutate(Total=.data$Effort*(lnorm.mean(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))-0.1),
                 TotalVar=.data$Effort^2*lnorm.se(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))^2)
        sim=replicate(nsim,rlnorm(nObs,
                                  meanlog=as.vector(a %*% mvrnorm(1,coefvals1,vcovvals1))+NewRandomVals,
                                  sdlog=sigma(modfit1))-0.1)*newdat$Effort
      }
      if(modtype  %in% c("Gamma","TMBgamma")) {
        if(class(modfit1)[1]=="glm") shapepar<-gamma.shape(modfit1)[[1]] else
          shapepar<-1/(glmmTMB::sigma(modfit1))^2
        allpred<-cbind(newdat,response1)   %>%
          mutate(Total=.data$Effort*(.data$fit-0.1),
                 TotalVar=.data$Effort^2*(.data$se.fit^2+.data$fit*shapepar))
        sim=replicate(nsim,newdat$Effort*(simulateGammaDraw(modfit1,nObs,a,NewRandomVals)-0.1) )
      }
      if(modtype %in% c("Delta-Lognormal","TMBdelta-Lognormal")) {
        allpred<-cbind(newdat,response1,response2) %>%
          mutate(pos.cpue=lnorm.mean(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)),
                 pos.cpue.se=lnorm.se(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)),
                 prob.se=sqrt(.data$se.fit^2+.data$fit*(1-.data$fit))) %>%
          mutate(Total=.data$Effort*.data$fit*.data$pos.cpue,
                 TotalVar=.data$Effort^2*lo.se(.data$fit,.data$prob.se,.data$pos.cpue,.data$pos.cpue.se)^2)
        sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coefvals1,vcovvals1) ) ))
        sim2=replicate(nsim,newdat$Effort*exp(rnorm(nObs,b %*%
                                                      mvrnorm(1,coefvals2,vcovvals2),sigma(modfit2)) ) )
        sim=sim1*sim2
      }
      if(modtype %in% c("Delta-Gamma","TMBdelta-Gamma")) {
        if(class(modfit2)[1]=="glm") shapepar<-gamma.shape(modfit2)[[1]] else
          shapepar<-1/(glmmTMB::sigma(modfit2))^2
        allpred<-cbind(newdat,response1,response2) %>%
          mutate(pos.cpue.se=sqrt(.data$se.fit2^2+.data$fit2*shapepar),
                 prob.se=sqrt(.data$se.fit^2+.data$fit*(1-.data$fit))) %>%
          mutate(Total=.data$Effort*.data$fit*.data$fit2,
                 TotalVar=.data$Effort^2*lo.se(.data$fit,.data$prob.se,.data$fit2,.data$pos.cpue.se)^2)
        sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coefvals1,vcovvals1)+NewRandomVals) ) )
        sim2=replicate(nsim,newdat$Effort*simulateGammaDraw(modfit2,nObs,b,NewRandomVals2) )
        sim=sim1*sim2
      }
      if(modtype=="NegBin") {
        allpred<-cbind(newdat,response1)   %>%
          mutate(Total=.data$fit,TotalVar=.data$se.fit^2+.data$fit+.data$fit^2/modfit1$theta)
        sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,coefvals1,vcovvals1))*newdat$Effort,
                                     size=modfit1$theta))  #Simulate negative binomial data
      }
      if(modtype=="Tweedie") {
        allpred=cbind(newdat,response1)   %>%
          mutate(Total=.data$Effort*.data$fit,
                 TotalVar=.data$Effort^2*(.data$se.fit^2+modfit1$phi*.data$fit^modfit1$p))
        sim=replicate(nsim,rtweedie(nObs,power=modfit1$p,
                                    mu=as.vector(exp(a %*% mvrnorm(1,coef(modfit1),modfit1$vcov))),
                                    phi=modfit1$phi))*newdat$Effort
      }
      if(modtype=="TMBnbinom1") {
        allpred<-cbind(newdat,response1)  %>%
          mutate(Total=.data$fit,
                 TotalVar=.data$se.fit^2+.data$fit+.data$fit*sigma(modfit1))
        sim = replicate(nsim,simulateNegBin1Draw(modfit1,nObs,a,newdat$Effort,NewRandomVals))
      }
      if(modtype=="TMBnbinom2") {
        allpred<-cbind(newdat,response1)  %>%
          mutate(Total=.data$fit,
                 TotalVar=.data$se.fit^2+.data$fit+.data$fit^2/sigma(modfit1))
        sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,fixef(modfit1)[[1]],
                                                               vcov(modfit1)[[1]])+NewRandomVals)*newdat$Effort, size=sigma(modfit1)))
      }
      if(modtype=="TMBtweedie") {
        allpred<-cbind(newdat,response1)  %>%
          mutate(Total=.data$Effort*.data$fit,
                 TotalVar=.data$Effort^2*(.data$se.fit^2+sigma(modfit1)*.data$fit^(glmmTMB::family_params(modfit1))))
        sim=replicate(nsim, simulateTMBTweedieDraw(modfit1,nObs,a,newdat$Effort,NewRandomVals) )
      }
      if(includeObsCatch & !modtype %in% c("Binomial","TMBbinomial")) {
        obsdatvalyear=obsdatval[obsdatval$Year==years[i],]
        d=match(allpred$matchColumn,obsdatvalyear$matchColumn)
        allpred$Total[!is.na(d)]= allpred$Total[!is.na(d)] + obsdatvalyear$Catch[d[!is.na(d)]]
        sim[!is.na(d),]= sim[!is.na(d),] + obsdatvalyear$Catch[d[!is.na(d)]]
      }
      if(includeObsCatch & modtype %in% c("Binomial","TMBbinomial")) {
        obsdatvalyear=obsdatval[obsdatval$Year==years[i],]
        d=match(allpred$matchColumn,obsdatvalyear$matchColumn)
        allpred$Total[!is.na(d)]= obsdatvalyear$pres[d[!is.na(d)]]
        sim[!is.na(d),]= obsdatvalyear$pres[d[!is.na(d)]]
      }
      stratatotal<-allpred %>%
        group_by_at(all_of(requiredVarNames)) %>%
        summarize(Total=sum(.data$Total,na.rm=TRUE))
      yeartotal<-allpred%>% group_by(.data$Year) %>%
        summarize(Total=sum(.data$Total,na.rm=TRUE))
      stratapredyear<-cbind(newdat,sim) %>%
        group_by_at(all_of(c("strata",requiredVarNames))) %>%
        summarize_at(.vars=as.character(1:nsim),.funs=sum,na.rm=TRUE) %>%
        rowwise() %>%
        mutate(Total.mean=mean(c_across(as.character(1:nsim))),
               TotalVar=var(c_across(as.character(1:nsim))),
               TotalLCI=quantile(c_across(as.character(1:nsim)),p=CIval/2),
               TotalUCI=quantile(c_across(as.character(1:nsim)),p=1-CIval/2)) %>%
        mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI),Total.mean=ifelse(.data$TotalLCI<0,0,.data$Total.mean)) %>%
        mutate(Total.se=sqrt(.data$TotalVar))  %>%
        mutate(Total.cv=.data$Total.se/.data$Total.mean)  %>%
        dplyr::select(-one_of(as.character(1:nsim)))
      stratapredyear$Total=stratatotal$Total
      yearpredyear<-cbind(newdat,sim) %>%
        group_by(.data$Year) %>%
        summarize_at(.vars=as.character(1:nsim),.funs=sum,na.rm=TRUE) %>%
        rowwise() %>%
        mutate(Total.mean=mean(c_across(as.character(1:nsim))),
               TotalVar=var(c_across(as.character(1:nsim))),
               TotalLCI=quantile(c_across(as.character(1:nsim)),p=CIval/2),
               TotalUCI=quantile(c_across(as.character(1:nsim)),p=1-CIval/2)) %>%
        mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI),Total.mean=ifelse(.data$Total.mean<0,0,.data$Total.mean)) %>%
        mutate(Total.se=sqrt(.data$TotalVar))  %>%
        mutate(Total.cv=.data$Total.se/.data$Total.mean)  %>%
        dplyr::select(-one_of(as.character(1:nsim)))
      yearpredyear$Total<-yeartotal$Total
      yearpred[i,c("Total.mean", "TotalVar", "TotalLCI", "TotalUCI", "Total.se" ,"Total.cv", "Total")]<-
        yearpredyear[1,c("Total.mean", "TotalVar", "TotalLCI", "TotalUCI", "Total.se" ,"Total.cv", "Total")]
      if(nrow(stratapredyear)>1) #If there are more than one strata in a year
        stratapred[stratapred$Year==years[i],c("strata","Total.mean", "TotalVar", "TotalLCI", "TotalUCI", "Total.se" ,"Total.cv", "Total")]<-
        stratapredyear[,c("strata","Total.mean", "TotalVar", "TotalLCI", "TotalUCI", "Total.se" ,"Total.cv", "Total")]
      if(nrow(stratapredyear)==1) #If there may be more than one year in a stratum (e.g. multiyear data)
        stratapred[stratapred$strata==stratapredyear$strata[1],c("Total.mean", "TotalVar","Total")]<-
        stratapred[stratapred$strata==stratapredyear$strata[1],c("Total.mean", "TotalVar","Total")] +
        stratapredyear[,c("Total.mean", "TotalVar", "Total")]
    }
  }
  if(is.na(max(yearpred$Total.cv)) | max(yearpred$Total.cv,na.rm=TRUE)>10) {
    print(paste(common[run],modtype," CV >10 or NA variance"))
    returnval=NULL
  }  else  {     returnval=yearpred  }
  if(printOutput) {
    write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modelScenario,modtype,"StratumSummary.csv"), row.names = FALSE)
    write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modelScenario,modtype,"AnnualSummary.csv"), row.names = FALSE)
  }
  returnval
}

#' Generate standard errors and confidence intervals of predictions with delta-method separately by year
#'
#' @param modfit1 Value
#' @param newdat Value
#' @param modtype Value
#' @param obsdatval Value
#' @param includeObsCatch Value
#' @param requiredVarNames Value
#' @param CIval Value
#' @param printOutput Value
#' @param catchType Value
#' @param common Value
#' @param dirname Value
#' @param run Value
#' @param modelScenario Value
#' @importFrom stats delete.response terms qnorm
#' @keywords internal
makePredictionsDeltaVar<-function(modfit1, newdat, modtype,  obsdatval, includeObsCatch,
                                  requiredVarNames, CIval, printOutput=TRUE, catchType, common, dirname, run,
                                  modelScenario) {
  if(modtype %in% c("Delta-Lognormal","Delta-Gamma","Tweedie","TMBdelta-Lognormal","TMBdelta-Gamma")) stop("No delta-method variance available, use simulute")
  #Separate out sample units
  if(includeObsCatch)    newdat$Effort=newdat$unsampledEffort/newdat$SampleUnits else
    newdat$Effort=newdat$Effort/newdat$SampleUnits
  newdat=uncount(newdat,.data$SampleUnits)
  nObs=nrow(newdat)
  newdat$SampleUnits=rep(1,nObs)
  newdatall=newdat
  #Set up output dataframes
  years=sort(unique(newdat$Year))
  yearpred=expand.grid(Year=years,Total=NA,TotalVar=NA)
  stratapred=data.frame(newdatall[!duplicated(newdatall$strata),requiredVarNames])
  stratapred$Total=stratapred$TotalVar=NA
  for(i in 1:length(years)) {
    newdat = newdatall[newdatall$Year==years[i],]
    #Get model matrix
    tm = delete.response(terms(modfit1))
    a = model.matrix(tm, newdat)
    ## predicted value
    if(modtype=="Tweedie" ) {
     predvallink = cplm::predict(modfit1,newdat=newdat)
     predval = cplm::predict(modfit1,newdat=newdat,type="response")

    } else {
    if(grepl("TMB",modtype)) {
      predvallink = predict(modfit1,newdat=newdat,allow.new.levels = TRUE)
      predval = predict(modfit1,newdat=newdat,type="response",allow.new.levels = TRUE)
    } else {
      predvallink = predict(modfit1,newdat=newdat)
      predval = predict(modfit1,newdat=newdat,type="response")
    }
    }
    if(grepl("TMB",modtype))  vcovval = a %*% vcov(modfit1)[[1]] %*% t(a) else
        vcovval = a %*% vcov(modfit1) %*% t(a)
    if(modtype %in% c("Binomial","TMBbinomial")) {
      residvar =  predval * (1-predval)  #Binomial variance
      deriv =  as.vector(exp(predvallink)/(exp(predvallink)+1)^2)
    }
    if(modtype %in% c("Normal","TMBnormal")) {
      predval=predval*newdat$Effort
      residvar =  rep(sigma(modfit1)^2,nrow(newdat))*newdat$Effort^2
      deriv =  rep(1,nrow(newdat))
    }
    if(modtype %in% c("Lognormal","TMBlognormal" )) {
      temp = predict(modfit1,newdata=newdat,se.fit=TRUE)
      predval = (lnorm.mean(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2))-0.1)*newdat$Effort
      # residvar = lnorm.se(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2))^2*newdat$Effort^2
      # deriv =  exp(temp$fit)*newdat$Effort
      deriv = (lnorm.mean(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2)))*newdat$Effort
      residvar = lnorm.se(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2))^2*newdat$Effort^2
    }
    if(modtype %in% c("Gamma","TMBgamma") ) {
      if(class(modfit1)[1]=="glm") shapepar<-gamma.shape(modfit1)[[1]] else
        shapepar<-1/(glmmTMB::sigma(modfit1))^2
      predval = (predval-0.1) * newdat$Effort
      predval[predval<0]<-0
      residvar = exp(predvallink)*newdat$Effort * shapepar
      deriv =  exp(predvallink) #derivative of exp(x) is exp(x)
    }
    if(modtype == c("NegBin") ) {
      residvar =  predval+predval^2/modfit1$theta
      deriv =  predval  #derivative of exp(x) is exp(x)
    }
    if(modtype %in% c("TMBnbinom1") ) {
      residvar =  predval+predval*sigma(modfit1)
      deriv =  predval  #derivative of exp(x) is exp(x)
    }
    if(modtype %in% c("TMBnbinom2") ) {
      residvar =  predval+predval^2/sigma(modfit1)
      deriv =  predval  #derivative of exp(x) is exp(x)
    }
    if(modtype=="Tweedie") {
      residvar =  (modfit1$phi*predval^modfit1$p)*newdat$Effort^2
      predval = predval * newdat$Effort
      deriv =  predval  #derivative of exp(x) is exp(x)
    }
    if(modtype=="TMBtweedie") {
      residvar =  sigma(modfit1)*predval^(glmmTMB::family_params(modfit1))*newdat$Effort^2
      predval = predval * newdat$Effort
      deriv =  predval  #derivative of exp(x) is exp(x)
    }
    if(includeObsCatch & ! modtype %in% c("Binomial","TMBbinomial")) {
      obsdatvalyear=obsdatval[obsdatval$Year==years[i],]
      d=match(newdatall$matchColumn,obsdatvalyear$matchColumn)
      #d=match(logdat$matchColumn,obsdatvalyear$matchColumn)
      d=d[!is.na(d)]
      predval[d]= predval[d] + obsdatvalyear$Catch
    }
    if(includeObsCatch & modtype %in% c("Binomial","TMBbinomial")) {
      obsdatvalyear=obsdatval[obsdatval$Year==years[i],]
      d=match(newdatall$matchColumn,obsdatvalyear$matchColumn)
      #d=match(logdat$matchColumn,obsdatvalyear$matchColumn)
      d=d[!is.na(d)]
      predval[d]= obsdatvalyear$pres
    }
    yearpred$Total[i]<-sum(predval)
    yearpred$TotalVar[i] = t(deriv) %*%
      vcovval %*%  deriv + sum(residvar)
    if(length(requiredVarNames)>1 )  {
      strata=unique(newdat$strata)
      for(j in 1:length(strata)) {
        stratapred$Total[stratapred$Year==years[i]][j] = sum(predval[newdat$strata==strata[j]])
        stratapred$TotalVarl[stratapred$Year==years[i]][j] = t(deriv[newdat$strata==strata[j]]) %*%
          vcovval[newdat$strata==strata[j],newdat$strata==strata[j]] %*%
          deriv[newdat$strata==strata[j]] + sum(residvar[newdat$strata==strata[j]])
      }
    }
  }
  if(length(requiredVarNames)>1) {
    stratapred<-stratapred %>%
      mutate(Total.se=sqrt(.data$TotalVar)) %>%
      mutate(Total.cv=.data$Total.se/.data$Total,
             Total.mean=NA,
             TotalLCI=.data$Total-qnorm(1-CIval/2)*.data$Total.se,
             TotalUCI=.data$Total+qnorm(1-CIval/2)*.data$Total.se) %>%
      mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI))
    if(printOutput)  write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modelScenario,modtype,"StratumSummary.csv"), row.names = FALSE)

  }
  yearpred<-yearpred %>%
    mutate(Total.se=sqrt(.data$TotalVar)) %>%
    mutate(Total.cv=.data$Total.se/.data$Total,
           Total.mean=NA,
           TotalLCI=.data$Total-qnorm(1-CIval/2)*.data$Total.se,
           TotalUCI=.data$Total+qnorm(1-CIval/2)*.data$Total.se)%>%
    mutate(TotalLCI=ifelse(.data$TotalLCI<0,0,.data$TotalLCI))
  if(is.na(max(yearpred$Total.cv)) | max(yearpred$Total.cv,na.rm=TRUE)>10) {
    print(paste(common[run],modtype," CV >10 or NA variance"))
    returnval=NULL
  }  else  {     returnval=yearpred  }
  if(printOutput) { #remove Total.mean column from csv
    write.csv(yearpred[,!(names(yearpred) %in% c("Total.mean"))],
              paste0(dirname[[run]],common[run],catchType[run],modelScenario,
                     modtype,"AnnualSummary.csv"), row.names = FALSE)
  }
  returnval
}


#' Generate predictions without estimating variance
#'
#' @param modfit1 Value
#' @param modfit2 Value
#' @param newdat Value
#' @param modtype Value
#' @param obsdatval Value
#' @param includeObsCatch Value
#' @param requiredVarNames Value
#' @param printOutput Value
#' @param catchType Value
#' @param common Value
#' @param dirname Value
#' @param run Value
#' @param modelScenario Value
#' @importFrom stats delete.response terms qnorm
#' @keywords internal
makePredictionsNoVar<-function(modfit1, modfit2=NULL, modtype, newdat, obsdatval=NULL,
                               nsims, includeObsCatch, requiredVarNames, printOutput=TRUE,
                               catchType, common, dirname, run,modelScenario) {
  if(includeObsCatch)    newdat$Effort=newdat$unsampledEffort/newdat$SampleUnits else
    newdat$Effort=newdat$Effort/newdat$SampleUnits
  newdat=uncount(newdat,.data$SampleUnits)
  getse=ifelse(modtype %in% c("Lognormal","Delta-Lognormal", "TMBlognormal","TMBdelta-Lognormal"),TRUE,FALSE)
  nObs=dim(newdat)[1]
  if(!is.null(modfit1)) {
    if(modtype=="Tweedie")  response1<-data.frame(cplm::predict(modfit1,newdata=newdat,type="response"))
    if(grepl("TMB",modtype))  response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=getse,allow.new.levels=TRUE))
    if(!grepl("TMB",modtype) & !modtype=="Tweedie")  response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=getse))
    if(dim(response1)[2]==1)     names(response1)="fit"
    if(!is.null(modfit2))  {
      if(grepl("TMB",modtype))  response2<-data.frame(predict(modfit2,newdata=newdat,type="response",se.fit=getse,allow.new.levels=TRUE))
      if(!grepl("TMB",modtype) )  response2<-data.frame(predict(modfit2,newdata=newdat,type="response",se.fit=getse))
      if(dim(response2)[2]==1) names(response2)[1]<-"fit"
      names(response2)=paste0(names(response2),"2")
    }
    if(modtype %in% c("Binomial","TMBbinomial")) {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=.data$fit)
    }
    if(modtype %in% c("Normal","TMBnormal","Tweedie","TMBtweedie")) {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=.data$fit*.data$Effort)
    }
    if(modtype %in% c("Lognormal","TMBlognormal")) {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=(lnorm.mean(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))-0.1)*.data$Effort)
    }
    if(modtype %in% c("Gamma","TMBgamma")) {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=(.data$fit-0.1)*.data$Effort)
    }
    if(modtype %in% c("Delta-Lognormal","TMBdelta-Lognormal" )) {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue=lnorm.mean(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2))) %>%
        mutate(Total=.data$Effort*.data$fit*.data$pos.cpue)
    }
    if(modtype %in% c("Delta-Gamma","TMBdelta-Gamma")) {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(Total=.data$Effort*.data$fit*.data$fit2)
    }
    if(modtype %in% c("NegBin","TMBnbinom1","TMBnbinom2")) {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$fit)
    }
    if(includeObsCatch & !modtype %in% c("Binomial","TMBbinomial")) {
      a=match(allpred$matchColumn,obsdatval$matchColumn)
      allpred$Total[!is.na(a)]= allpred$Total[!is.na(a)] + obsdatval$Catch[a[!is.na(a)]]
    }
    if(includeObsCatch & modtype %in% c("Binomial","TMBbinomial")) {
      a=match(allpred$matchColumn,obsdatval$matchColumn)
      allpred$Total[!is.na(a)]= obsdatval$pres[a[!is.na(a)]]
    }
    stratapred<-allpred %>%
      group_by_at(all_of(requiredVarNames)) %>%
      summarize(Total=sum(.data$Total,na.rm=TRUE)) %>%
      mutate(Total.mean=NA,TotalVar=NA,	TotalLCI=NA,	TotalUCI=NA,	Total.se=NA,
             Total.cv=NA)

    yearpred<-allpred%>% group_by(.data$Year) %>%
      summarize(Total=sum(.data$Total,na.rm=TRUE)) %>%
      mutate(Total.mean=NA,TotalVar=NA,	TotalLCI=NA,	TotalUCI=NA,	Total.se=NA,
             Total.cv=NA)
    returnval=yearpred
    if(printOutput) {
      write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modelScenario,modtype,"StratumSummary.csv"), row.names = FALSE)
      write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modelScenario,modtype,"AnnualSummary.csv"), row.names = FALSE)
    }
  } else {
    returnval=NULL
  }
  returnval
}


#' Function to get an abundance index with SE
#'
#' @param modfit1 Value
#' @param modfit2 Value
#' @param newdat Value
#' @param modType Value
#' @param nsims Value
#' @param printOutput Value
#' @param catchType Value
#' @param common Value
#' @param dirname Value
#' @param run Value
#' @param modelScenario Value
#' @keywords internal
makeIndexVar<-function(modfit1, modfit2=NULL, modType, indexVarNames,newdat, nsims,
                       printOutput=FALSE, catchType = NULL, common = NULL, dirname = NULL, run = NULL,
                       modelScenario) {
  returnval=NULL
  if(!is.null(modfit1)) {
    if(modType=="Tweedie")    response1<-data.frame(cplm::predict(modfit1,newdata=newdat,type="response",se.fit=TRUE)) else
    if(!grepl("TMB",modType)) response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(grepl("TMB",modType))  response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE,allow.new.levels=TRUE))
    if(dim(response1)[2]==1) {
      names(response1)="fit"
      if(modType=="Tweedie")
        response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp", offsetval=NULL, nsim=nsims) else
          response1$se.fit=rep(NA,dim(response1)[2])
    }
    if(!is.null(modfit2))  {
      if(!grepl("TMB",modType)) response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
      if(grepl("TMB",modType)) response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response",allow.new.levels=TRUE))
      names(response2)=paste0(names(response2),"2")
    }
    if(modType %in% c("Delta-Lognormal","TMBdelta-Lognormal" )) {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue=lnorm.mean(.data$fit2,.data$se.fit2),
               pos.cpue.se=lnorm.se(.data$fit2,.data$se.fit2),
               prob.se=.data$se.fit) %>%
        mutate(Index=.data$fit*.data$pos.cpue,
               SE=lo.se(.data$fit,.data$prob.se,.data$pos.cpue,.data$pos.cpue.se))
    }
    if(modType %in% c("Delta-Gamma","TMBdelta-Gamma")){
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue.se=.data$se.fit2,
               prob.se=.data$se.fit) %>%
        mutate(Index=.data$fit*.data$fit2,
               SE=lo.se(.data$fit,.data$prob.se,.data$fit2,.data$pos.cpue.se))
    }
    if(modType %in% c("Binomial","NegBin","Tweedie","TMBnbinom1","TMBnbinom2","Normal","TMBtweedie","TMBbinomial","TMBnormal")) {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Index=.data$fit, SE=.data$se.fit)
    }
    if(modType %in% c("Gamma","TMBgamma")) {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Index=.data$fit-0.1, SE=.data$se.fit)
    }
    if(modType %in% c("Lognormal","TMBlognormal")) {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Index=lnorm.mean(.data$fit,.data$se.fit)-0.1,
               SE=lnorm.se(.data$fit,.data$se.fit))
    }
    allpred=allpred %>%
      mutate(lowerCI=.data$Index-.data$SE,upperCI=.data$Index+.data$SE)  %>%
      mutate(lowerCI=ifelse(.data$lowerCI<0,0,.data$lowerCI))
    returnval=allpred
    if(printOutput) {
      indexVarNames <- indexVarNames[indexVarNames %in% colnames(allpred)]
      colsToSave <- c(indexVarNames, "Index","SE","lowerCI","upperCI")
      write.csv(allpred[, colsToSave, drop = FALSE],
                paste0(dirname[[run]],common[run],catchType[run],modelScenario,modType,"Index.csv"), row.names = FALSE)
    }
  }
  returnval
}

#' Function plots residuals with both R and Dharma library and calculate residual diagnostics.
#'
#' @param modfit1 Value
#' @param modType Value
#' @param fileName Value
#' @param nsim Value
#' @import quantreg grDevices DHARMa
#' @importFrom stats residuals qunif
#' @importFrom gridExtra grid.arrange
#' @keywords internal
ResidualsFunc<-function(modfit1,modType,fileName=NULL,nsim=250) {
  if(!is.null(fileName))   pdf(fileName,height=5,width=7)
  if(!is.null(modfit1)) {
    if(modType=="Tweedie")  dfcheck<-data.frame(Expected=cplm::predict(modfit1),Residuals=residuals(modfit1)) else
     dfcheck<-data.frame(Expected=predict(modfit1),Residuals=residuals(modfit1))
    if(nrow(dfcheck)>1000)
      subsam<-sample(1:nrow(dfcheck),1000) else
        subsam<-1:nrow(dfcheck) #Only show 1000 residuals if have more than that
      g1<-ggplot(dfcheck[subsam,],aes(x=.data$Expected,y=.data$Residuals))+geom_point()+
        geom_abline(intercept=0,slope=0)+
        ggtitle(paste("a. ",modType,"ordinary residuals"))
      g2<-ggplot(dfcheck[subsam,],aes(sample=.data$Residuals))+geom_qq()+geom_qq_line()+
        ggtitle(paste("b. QQ normal of residuals"))
      if(class(modfit1)[1] =="cpglm") {  #Extra step to simulate DHARMa for cpglm or mgcv
        simvals=simulateTweedie(modfit1,nsim)
        simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                            fittedPredictedResponse = cplm::predict(modfit1,type="response")))
        if(class(simulationOutput)[1]=="try-error") {
          simvals=simulateTweedie(modfit1,nsim*4)
          simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                              fittedPredictedResponse = cplm::predict(modfit1,type="response")))
        }
        if(class(simulationOutput)[1]=="try-error") {
          simvals=simulateTweedie(modfit1,nsim*10)
          simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                              fittedPredictedResponse = cplm::predict(modfit1,type="response")))
        }
      }
      if(class(modfit1)[1]=="gam" & modType=="NegBin") {
        simvals=simulateNegBinGam(modfit1,nsim)
        simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                            fittedPredictedResponse = predict(modfit1,type="response")))
        if(class(simulationOutput)[1]=="try-error") {
          simvals=simulateNegBinGam(modfit1,nsim*4)
          simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                              fittedPredictedResponse = predict(modfit1,type="response")))
        }
        if(class(simulationOutput)[1]=="try-error") {
          simvals=simulateNegBinGam(modfit1,nsim*10)
          simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                              fittedPredictedResponse = predict(modfit1,type="response")))
        }
        if(class(simulationOutput)[1]=="try-error") {
          simvals=simulateNegBinGam(modfit1,nsim*15)
          simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                              fittedPredictedResponse = predict(modfit1,type="response")))
        }
      }
      if( class(modfit1)[1] !="cpglm" & !(class(modfit1)[1]=="gam" & modType=="NegBin"))     {  #Regular DHARMa residuals work for everything but cpglm
        simulationOutput <- try(simulateResiduals(fittedModel = modfit1, n = nsim))
        if(class(simulationOutput)[1]=="try-error") simulationOutput <- try(simulateResiduals(fittedModel = modfit1, n = nsim*4))
        if(class(simulationOutput)[1]=="try-error") simulationOutput <- try(simulateResiduals(fittedModel = modfit1, n = nsim*10))
      }
      if(class(simulationOutput)[1]!="try-error") {
        #      plot(simulationOutput, quantreg = F)
        #      title(modType,outer=2,line=-1)
        df1<-data.frame(Residual=simulationOutput$scaledResiduals[subsam],Predictor=simulationOutput$fittedPredictedResponse[subsam]) %>%
          arrange(.data$Residual) %>%
          mutate(Empirical.Quantile=(1:n())/(n()+1)) %>%
          mutate(Expected.Quantile=qunif(.data$Empirical.Quantile)) %>%
          mutate(Rank.Predictor= rank(.data$Predictor, ties.method = "average")) %>%
          mutate(Rank.Predictor = .data$Rank.Predictor/max(.data$Rank.Predictor))
        g3<-ggplot(df1,aes(x=.data$Expected.Quantile,y=.data$Residual))+geom_point()+
          geom_abline()+ylab("DHARMa scaled residuals")+ xlab("Expected quantile")+
          ggtitle("c. QQ uniform scaled residuals")
        if(length(unique(df1$Rank.Predictor))>1) {
          g4<-ggplot(df1,aes(x=.data$Rank.Predictor,y=.data$Residual))+
            geom_point()+xlab("Model predictions (rank transformed)")+
            ylab("DHARMa scaled residuals")+ggtitle("d. Scaled residual vs. predicted")+
            geom_hline(aes(yintercept=0.5),lty=2)+
            geom_hline(aes(yintercept=0.75),lty=2)+
            geom_hline(aes(yintercept=0.25),lty=2)
          if(class(try(rqss(Residual~qss(Rank.Predictor,lambda=2),data=df1),silent = TRUE))!="try-error")
            g4<-g4+geom_quantile(method = "rqss",col="red", formula=y ~ qss(x, lambda = 2))
        } else {
          g4<-ggplot(df1,aes(x=.data$Rank.Predictor,y=.data$Residual))+
            geom_point()+xlab("Model predictions (rank transformed)")+
            ylab("DHARMa scaled residuals")+ggtitle("d. Scaled residual vs. predicted")+
            geom_hline(aes(yintercept=0.5),lty=2)+geom_hline(aes(yintercept=0.75),lty=2)+geom_hline(aes(yintercept=0.25),lty=2)
        }
        grid.arrange(g1,g2,g3,g4,ncol=2)
        test1=testUniformity(simulationOutput,plot=FALSE)
        test2=testDispersion(simulationOutput,plot=FALSE)
        test3=testZeroInflation(simulationOutput,plot=FALSE)
        test4=suppressWarnings(testOutliers(simulationOutput,plot=FALSE))
        returnval=c(test1$statistic,
                    test1$p.value,
                    test2$statistic,
                    test2$p.value,
                    test3$statistic,
                    test3$p.value,
                    test4$statistic,
                    test4$p.value)
        names(returnval)=c("KS.D","KS.p","Dispersion.ratio","Dispersion.p","ZeroInf.ratio",
                           "ZeroInf.p","Outlier","Outlier.p")
        if(modType %in% c("Delta-Gamma","Delta-Lognormal","Lognormal")) returnval[5:6]=NA
      } else returnval=NULL
  } else returnval=NULL
  if(!is.null(fileName))  dev.off()
  returnval
}

#' Function to fit a specified model formula and print outputs for Cross validation
#'
#' @param formula1 Value
#' @param modType Value
#' @param obsdatval Value
#' @keywords internal
FitModelFuncCV<-function(formula1,modType,obsdatval) {
  if(modType %in% c("Binomial") )  {
    obsdatval$y=obsdatval$pres
    modfit1<-try(glm(formula1,data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=1000)))
  }
  if(modType %in% c("TMBbinomial") )  {
    obsdatval$y=obsdatval$pres
    modfit1<-try(glmmTMB(formula1,data=obsdatval,family="binomial"))
  }
  if(modType=="Normal") {
    obsdatval$y=obsdatval$cpue
    modfit1<-try(lm(formula1,data=obsdatval))
  }
  if(modType=="Gamma") {
    obsdatval$y=obsdatval$cpue+0.1
    modfit1<-try(glm(formula1,data=obsdatval,family=Gamma(link=log)))
  }
  if(modType=="TMBgamma") {
    obsdatval$y=obsdatval$cpue+0.1
    modfit1<-try(glmmTMB(formula1,data=obsdatval,family=Gamma(link=log)))
  }
  if(modType=="TMBnormal") {
    obsdatval$y=obsdatval$cpue
    modfit1<-try(glmmTMB(formula1,data=obsdatval))
  }
  if(modType=="Lognormal") {
    obsdatval$y=log(obsdatval$cpue+0.1)
    modfit1<-try(lm(formula1,data=obsdatval))
  }
  if(modType=="TMBlognormal") {
    obsdatval$y=log(obsdatval$cpue+0.1)
    modfit1<-try(glmmTMB(formula1,data=obsdatval))
  }
  if(modType=="Delta-Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(lm(formula1,data=obsdatval))
  }
  if(modType=="TMBdelta-Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(glmmTMB(formula1,data=obsdatval))
  }
  if(modType=="Delta-Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(glm(formula1,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="TMBdelta-Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(glmmTMB(formula1,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="NegBin") {
    obsdatval$y=round(obsdatval$Catch)
    modfit1=try(glm.nb(formula1,data=obsdatval,control=glm.control(epsilon=1E-6,maxit=45),na.action=na.fail))
  }
  if(modType=="Tweedie") {
    obsdatval$y=obsdatval$cpue
    modfit1=try(cplm::cpglm(formula1,data=obsdatval))
  }
  if(modType %in% c("TMBnbinom1","TMBnbinom2") ){
    obsdatval$y=round(obsdatval$Catch)
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula1,family=TMBfamily,data=obsdatval))
  }
  if(modType =="TMBtweedie"){
    obsdatval$y=obsdatval$cpue
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula1,family=TMBfamily,data=obsdatval))
  }
  if(class(modfit1)[1]=="try-error") modfit1=NULL
  modfit1
}

#' Function to predict CPUE without variances to get predictions quickly for cross validation
#'
#' @param modfit1 Value
#' @param modfit2 Value
#' @param modType Value
#' @param newdat Value
#' @param obsdatval Value
#' @keywords internal
makePredictions<-function(modfit1,modfit2=NULL,modType,newdat,obsdatval=NULL) {
  if(!is.null(modfit1)) {
    if(modType=="Tweedie")    predval1<-try(data.frame(fit=cplm::predict(modfit1,newdata=newdat,type="response"))) else
      predval1<-try(data.frame(predict(modfit1,newdata=newdat,se.fit=TRUE,type="response")))
    if(class(predval1)[[1]]!="try-error") {
      if(!is.null(modfit2))  {
        predval2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
        names(predval2)=paste0(names(predval2),"2")
      }
      if(modType %in% c("Delta-Lognormal","TMBdelta-Lognormal")) {
        allpred<-cbind(newdat,predval1,predval2)  %>%
          mutate(est.cpue=.data$fit*lnorm.mean(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)))
      }
      if(modType %in% c("Delta-Gamma","TMBdelta-Gamma")) {
        allpred<-cbind(newdat,predval1,predval2)   %>%
          mutate(est.cpue=.data$fit*.data$fit2)
      }
      if(modType %in% c("NegBin","TMBnbinom1","TMBnbinom2")) {
        allpred<-cbind(newdat,predval1)   %>%
          mutate(est.cpue=.data$fit/.data$Effort)
      }
      if(modType %in% c("TMBtweedie","Normal","TMBnormal","Tweedie")){
        allpred<-cbind(newdat,predval1)   %>%
          mutate(est.cpue=.data$fit)
      }
      if(modType %in% c("Gamma","TMBgamma")) {
        allpred<-cbind(newdat,predval1)   %>%
          mutate(est.cpue=.data$fit-0.1)
      }
      if(modType %in% c("Lognormal","TMBlognormal")){
        allpred<-cbind(newdat,predval1)   %>%
          mutate(est.cpue=lnorm.mean(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))-0.1)
      }
      returnval=allpred
    } else returnval=NULL
  } else {
    returnval=NULL
  }
  returnval
}


#' Function to simulate DHARMa residuals from a negative binomial GAM or GLM
#'
#' @param modfit Value
#' @param nsims Value
#' @param offsetval Value
#' @importFrom stats residuals predict
#' @keywords internal
simulateNegBinGam <- function(modfit, nsims=250, offsetval=1){
  muval = predict(modfit, type = "response")*offsetval  #Get the mean with offset
  nObs = length(muval)
  thetaval = modfit$family$getTheta(trans=TRUE)  #Get theta not log transforme
  sim = replicate(nsims,rnbinom(nObs,mu=muval, size=thetaval))  #Simulate negative binomial data
  sim
}


#' Generate standard errors of predictions from simulation from regression coefficients and their var/covar matrix
#'
#' @param modfit Value
#' @param df1 Value
#' @param transFunc Value
#' @param offsetval Value
#' @param nsim Value
#' @import tidyr
#' @importFrom stats predict model.matrix
#' @importFrom MASS mvrnorm
#' @keywords internal
getSimSE<-function(modfit, df1, transFunc="none", offsetval=NULL, nsim) {
  if(class(modfit)[1]=="cpglm") vcovval=modfit$vcov else
      vcovval=vcov(modfit)
  if(length(coef(modfit))==dim(vcovval)[1]) {
    b=t(mvrnorm(nsim,coef(modfit),vcovval))
    yvar=sub( " ", " ",formula(modfit) )[2]
    df1<-cbind(y=rep(1,dim(df1)[1]),df1)
    names(df1)[1]<-yvar
    a=model.matrix(formula(modfit),data=df1)
    d<- a %*% b
    if(transFunc=="exp") d<-exp(d)
    if(transFunc=="ilogit") d<-ilogit(d)
    e<-apply(d,1,sd)
    if(!is.null(offsetval)) e<-e*data.frame(df1)[,offsetval]  #For negbin only
  } else {
    e<-rep(NA,dim(df1)[1])
  }
  e
}

#' Inverse logit
#' @keywords internal
ilogit=function(x) {
  1/(1+exp(-x))
}

#' Calculate lognormal mean and standard error from normal mean and se
#' @keywords internal
lnorm.mean=function(x1,x1e) {
  exp(x1+0.5*x1e^2)
}

#' lnorm.se
#' @keywords internal
lnorm.se=function(x1,x1e) {
  ((exp(x1e^2)-1)*exp(2*x1+x1e^2))^0.5
}

#' as.numeric2
#' @keywords internal
as.numeric2<-function(x) {
  as.numeric(as.character(x))
}

#' simulateNegBin1Draw
#'
#' @param modfit Value
#' @param nObs Value
#' @param b Value
#' @param Effort Value
#' @param NewRandomVals Value
#' @importFrom MASS mvrnorm
#' @importFrom glmmTMB fixef
#' @importFrom stats sigma rnbinom vcov
#' @keywords internal
simulateNegBin1Draw<-function(modfit,nObs,b,Effort,NewRandomVals=0) {
  muval<-exp(b %*% mvrnorm(1,fixef(modfit)[[1]],vcov(modfit)[[1]]) + NewRandomVals)*Effort
  thetaval<-sigma(modfit)
  predval<-rep(0,nObs)
  predval[muval>0]<-rnbinom(length(muval[muval>0]),mu=muval[muval>0],size=muval[muval>0]/thetaval)
  predval
}

#' simulateGammaDraw
#'
#' @param modfit Value
#' @param nObs Value
#' @param b Value
#' @param NewRandomVals Value
#' @importFrom MASS mvrnorm gamma.shape
#' @importFrom glmmTMB fixef
#' @importFrom stats coef vcov rgamma
#' @keywords internal
simulateGammaDraw<-function(modfit,nObs,b,NewRandomVals) {
  if(class(modfit)[1]=="glm") {
   muval<-exp(b %*% mvrnorm(1,coef(modfit),vcov(modfit))+NewRandomVals)
   shapeval<-gamma.shape(modfit)[[1]]
  } else {
    muval<-exp(b %*% mvrnorm(1,fixef(modfit)[[1]],vcov(modfit)[[1]]))
    shapeval<-1/glmmTMB::sigma(modfit)^2
  }
  scaleval<-muval/shapeval
  rgamma(nObs,shape=shapeval,scale=scaleval)
}

#' simulateTMBTweedieDraw
#'
#' @param modfit Value
#' @param nObs Value
#' @param b Value
#' @param Effort Value
#' @param NewRandomVals Value
#' @importFrom MASS mvrnorm
#' @importFrom stats vcov rgamma sigma
#' @importFrom tweedie rtweedie
#' @keywords internal
simulateTMBTweedieDraw<-function(modfit,nObs,b,Effort,NewRandomVals=0) {
  muval<-as.vector(exp(b %*% mvrnorm(1,fixef(modfit)[[1]],vcov(modfit)[[1]])+NewRandomVals))
  if(all(muval>0)) {
    simval<-rtweedie(nObs,power=glmmTMB::family_params(modfit),
                     mu=muval,phi=sigma(modfit))*Effort  } else  {
                       simval=rep(NA,nObs)
                     }
  simval
}

#' Generate simulations to use as input to DHARMa residual calculations for the Tweedie from cpglm.
#'
#' @param modfit1 Value
#' @param nsims Value
#' @importFrom stats predict
#' @importFrom tweedie rtweedie
#' @keywords internal
simulateTweedie <- function(modfit1, nsims){
  pred = cplm::predict(modfit1, type = "response")
  nObs = length(pred)
  sim = replicate(nsims,rtweedie(nObs,xi=modfit1$p, mu=pred,phi=modfit1$phi))
  return(sim)
}

#' Generate mean and standard error of predictions for delta lognormal by simulation
#'
#' @param modfitBin Value
#' @param modfitLnorm Value
#' @param df1 Value
#' @param nsim Value
#' @importFrom stats coef vcov model.matrix
#' @importFrom MASS mvrnorm
#' @keywords internal
getSimDeltaLN<-function(modfitBin,modfitLnorm, df1, nsim=10000) {
  b1=t(mvrnorm(nsim,coef(modfitBin),vcov(modfitBin)))
  yvar=sub( " ", " ",formula(modfitBin) )[2]
  df11<-cbind(y=rep(1,dim(df1)[1]),df1)
  names(df11)[1]<-yvar
  a1=model.matrix(formula(modfitBin),data=df11)
  d1<- a1 %*% b1
  b2=t(MASS::mvrnorm(nsim,coef(modfitLnorm),vcov(modfitLnorm)))
  yvar=sub( " ", " ",formula(modfitLnorm) )[2]
  df12<-cbind(y=rep(1,dim(df1)[1]),df1)
  names(df12)[1]<-yvar
  a2=model.matrix(formula(modfitLnorm),data=df12)
  d2<- a2 %*% b2
  lnormmean<-apply(exp(d2),1,mean)
  lnormse<-apply(exp(d2),1,sd)
  deltamean<-apply(ilogit(d1)*exp(d2),1,mean)
  deltase<-apply(ilogit(d1)*exp(d2),1,sd)
  data.frame(delta.mean=deltamean,delta.se=deltase,lnorm.mean=lnormmean,lnorm.se=lnormse)
}


#######Delta lognormal functions

#' Variance of a product from Lo et al. (1992)
#'
#' @param x1 Value
#' @param x1e Value
#' @param x2 Value
#' @param x2e Value
#' @importFrom stats cor
#' @keywords internal
lo.se=function(x1,x1e,x2,x2e) {
  if(length(x1[!is.na(x1) &!is.na(x2)])>1 & length(unique(x1))>1 & length(unique(x2))>1 )
   cor12=cor(x1[!is.na(x1) &!is.na(x2)],x2[!is.na(x1) &!is.na(x2)]) else cor12=0
  (x1e^2 * x2^2 +x2e^2*x1^2+2*x1*x2*cor12*x1e*x2e)^0.5
}

#' Calculate normal mean from lognormal mean and se
#'
#' @param x1 Value
#' @param x1e Value
#' @keywords internal
norm.mean=function(x1,x1e) {
  log(x1)-0.5*log(1+(x1e/x1)^2)
}

#' Calculate standard error from lognormal mean and se
#'
#' @param x1 Value
#' @param x1e Value
#' @keywords internal
norm.se=function(x1,x1e) {
  sqrt(log(1+(x1e/x1)^2))
}

#' Calculate RMSE
#' @keywords internal
getRMSE<-function(yhat,y) {
  if(!is.null(yhat))
    rmse=sqrt(sum((yhat-y)^2)/length(y)) else
      rmse=NA
    rmse
}

#' Calculate mean error
#' @keywords internal
getME<-function(yhat,y) {
  if(!is.null(yhat))
    me=sum(yhat-y)/length(y) else
      me=NA
    me
}

#----------------------------


# Statified mean estimator (Cochran)
# stratified.func=function(y,g,N) {
#   ymean=tapply(y,g,mean)
#   n=tapply(x,g,length)
#   yvar=tapply(y,g,var)
#   stratum.est=ymean*N
#   stratum.var=N*(N-n)*yvar/n
#   total.est=sum(stratum.est)
#   total.var=sum(stratum.var)
#   list(stratum.est=stratum.est,stratum.se=sqrt(stratum.var),total.est=sum(stratum.est))
# }

#' Function to convert data in excel format with date and time separated by a blank into an R format date
#'
#' @param x data
#' @param dateformat Date format, default is %m/%d/%Y.
#' @importFrom reshape2 colsplit
#' @keywords internal
getdatefunc=function(x, dateformat="%m/%d/%Y") {
  y=reshape2::colsplit(as.character(x)," ",names=c("date","time"))
  as.Date(y[,1],format=dateformat)
}


#' Function 2 to convert data in excel format with date and time separated by a blank into an R format date
#'
#' @param x data
#' @param dateformat Date format, default is %d%b%Y.
#' @importFrom reshape2 colsplit
#' @keywords internal
getdatefunc2=function(x,dateformat="%d%b%Y") {
  y=reshape2::colsplit(as.character(x),":",names=c("date","hour","sec"))
  as.Date(y[,1],format=dateformat)
}


#' Convert time in excel format
#'
#' @param x data
#' @importFrom reshape2 colsplit
#' @keywords internal
gettimefunc2=function(x) {
  y=reshape2::colsplit(as.character(x),":",names=c("date","hour","min","sec"))
  timeval=y[,2]+y[,3]/60+y[,4]/3600
  timeval
}


#' Function to convert months into 2, 3 4 or 6 numbered seasons
#' @keywords internal
seasonfunc<-function(month,numseason=4) {
  if(!is.numeric(month)) month=as.numeric(as.character(month))
  seasons=rep(1:numseason,each=12/numseason)
  seasons[month]
}


#' Function to group the specified number of years together, starting from the last year
#' @keywords internal
yearfunc<-function(year,numyears=1) {
  if(!is.numeric(year)) year=as.numeric(as.character(year))
  trunc((year-max(year))/numyears)*numyears+max(year)
}


#' Function to look for positive and zero observations across levels of multiple factors
#' @keywords internal
CheckForPositives<-function(datval,species,variables) {
  tables=list()
  for(i in 1:length(variables)) {
    tables[[i]]=table(ifelse(datval[,species]>0,"Positive","Zero"),datval[,variables[i]])
    names(tables)[i]=variables[i]
  }
  tables
}

#' Same for plotting
#'
#' @param datval Value
#' @param species Value
#' @param variables Value
#' @import dplyr stringr
#' @keywords internal
CheckForPositivesPlot<-function(datval,species,variables) {
  tables=list()
  datval$Positive=ifelse(datval[,species]>0,"Positive","Zero")
  datval=select(datval,all_of(c("Positive",variables)))
  for(i in 1:length(variables)) {
    tables[[i]]=datval %>% rename(Group=!!variables[i]) %>%
      group_by(.data$Group,.data$Positive) %>%
      summarize(Count=n()) %>%
      mutate(Group=as.character(.data$Group))
  }
  names(tables)=variables
  tables=bind_rows(tables, .id="Variable")
  tables$Group=factor(tables$Group)
  tables$Group=factor(tables$Group,levels=str_sort(levels(tables$Group),numeric=TRUE))
  tables
}

#' Function to count the number of unique levels in a vector
#' @keywords internal
length.unique=function(x) length(unique(x))

#'Function to divide up areas. Input grid areas, returns East vs. West
#' @keywords internal
areaDivide<-function(area) {
  EW=ifelse(area>=11,"W","E")
  EW
}

#' Stratum designations from Scott-Denton paper, and shrimp observer manual for GOM shrimp areas
#' @keywords internal
areafunc=function(x) {
  x$StatZone[is.na(x$StatZone)]=-1
  area=rep(-1,dim(x)[1])
  area[x$StatZone>=1 & x$StatZone<=9]=1 #"WFL"
  area[x$StatZone>=10 & x$StatZone<=12]=2 #"AL-MS"
  area[x$StatZone>=13 & x$StatZone<=17]=3 #"LA"
  area[x$StatZone>=18 &x$StatZone<=21]=4 #"TX"
  x$LatInS[is.na(x$LatInS)]=0
  x$LatInM[is.na(x$LatInM)]=0
  x$LatOutS[is.na(x$LatOutS)]=0
  x$LatOutM[is.na(x$LatOutM)]=0
  x$lat=x$LatInD+x$LatInM/60+x$LatInS/3600
  x$lat[is.na(x$lat)]=(x$LatOutD+x$LatOutM/60+x$LatOutS/3600)[is.na(x$lat)]
  area[x$StatZone>=24 & x$StatZone<=29 ]=5 #"EFL"
  area[x$StatZone==30 & x$lat<30.708]=5 #"EFL"
  area[x$StatZone==30 & is.na(x$lat)]=5 #"EFL"
  area[x$StatZone>=30 & x$lat>=30.708 ]=6 #"GA"
  area[x$StatZone==31 ]=6 #"GA"
  area[x$StatZone==32 ]=7 #"SC"
  area[x$StatZone==33 & x$lat<33.86]=7 #"SC"
  area[x$StatZone==33 & x$lat>=33.86]=8 #"NC"
  area[x$StatZone>=34]=8 #"NC"
  area[is.na(x$StatZone) & x$LonInD>88 & x$LonInD<89]=2
  area
}


#' Function to find range of a numerical variable
#' @keywords internal
getRange<-function(x) {
  max(x,na.rm=TRUE)-min(x,na.rm=TRUE)
}

#' Function to convert new areas to old areas from Kevin McCarthy
#' @keywords internal
areaGOM=function(x) {
  x[x>=2383& x<= 2384]=2
  x[x >= 2483 & x <= 2485]=2
  x[x >= 2581 & x <= 2585]=3
  x[x >= 2681 & x <= 2685]=4
  x[x >= 2782 & x <= 2785]=5
  x[x >= 2882 & x <= 2884]=6
  x[x >= 2982 & x <= 2984]=7
  x[x >= 3083 & x <= 3084]=7
  x[x==3085 | x==2985 | x==2885]=8
  x[x %in% c(3086,2986,2886,2786,2686,2586,2486)] = 9
  x[x %in% c(3087,2987,2887,2787,2687,2587)] = 10
  x[x %in% c(3088,2988,2888,2788,2688,2588)] = 11
  x[x %in% c(3089,3090)] = 12
  x[x %in% c(2989,2889,2789,2689,2589)] = 13
  x[x %in% c(2990,2890,2790,2690,2590)] = 14
  x[x %in% c(2991,2891,2791,2691,2591)] = 15
  x[x %in% c(2992,2892,2792,2692,2592)] = 16
  x[x %in% c(2993,2893,2793,2693,2593)] = 17
  x[x %in% c(2994,2894,2794,2694,2594)] = 18
  x[x %in% c(2995,2895,2896)] = 19
  x[x %in% c(2797,2796,2795)] = 20
  x[x %in% c(2697,2696,2695)] = 21
  x[x>1000]=NA
  x
}


#' Function to calculate number of hours fished in each day counted as first set to last haul
#' @keywords internal
fishTimeFunc<-function(timeout,timein,prop.sampled=1) {
  timeout=ifelse(timeout>timein & !is.na(timein+timeout),timeout,timeout+24)
  maxtimeday=(max(timeout,na.rm=TRUE)-min(timein,na.rm=TRUE))
  if(maxtimeday<0.5) maxtimeday=0.5   #Arbitrary minimum
  meanprop=mean(prop.sampled,na.rm=TRUE)
  if(is.na(meanprop)) meanprop=1  #based on median proportion sampled=1
  maxtimeday*meanprop
}

#' fishTimeFunc
#' @keywords internal
fishTimeFunc<-function(timeout,timein,prop.sampled=1) {
  timeout=ifelse(timeout>timein & !is.na(timein+timeout),timeout,24)
  maxtimeday=(max(timeout,na.rm=TRUE)-min(timein,na.rm=TRUE))
  if(maxtimeday<0.5) maxtimeday=0.5   #Arbitrary minimum
  meanprop=mean(prop.sampled,na.rm=TRUE)
  if(is.na(meanprop)) meanprop=1  #based on median proportion sampled=1
  maxtimeday*meanprop
}

#' Exact variance of the product of two independent variables, from Goodman (1960)
#' @keywords internal
goodman.var<-function(x,y) {
  var(x)*y+var(y)*x-var(x)*var(y)
}

#' Function to plot either total positive trips (binomial) or total catch/bycatch (all other models)
#' @param yearpred Value
#' @param modType Value
#' @param fileName Value
#' @param subtext Value
#' @param allVarNames Value
#' @param startYear Value
#' @param common Value
#' @param run Value
#' @param catchType Value
#' @param catchUnit Value
#' @import dplyr
#' @keywords internal
plotSums<-function(yearpred,modType,fileName, subtext="", allVarNames, startYear, common, run, catchType, catchUnit,VarCalc) {
  if(is.numeric(yearpred$Year) & "Year" %in% allVarNames)
    yearpred$Year[yearpred$Year<startYear]=yearpred$Year[yearpred$Year<startYear]+startYear
 #  yearpred$Year[yearpred$Source!="Ratio"]=yearpred$Year[yearpred$Source!="Ratio"]+startYear
  if(!is.null(yearpred)) {
    if(modType=="Binomial") ytitle=paste0(common[run]," ","predicted total positive trips") else
      ytitle=paste0("Total ",common[run]," ",catchType[run]," (",catchUnit[run],")")
    yearpred<-yearpred %>%
      mutate(Year=as.numeric(as.character(.data$Year)),ymin=.data$Total-.data$Total.se,ymax=.data$Total+.data$Total.se) %>%
      mutate(ymin=ifelse(.data$ymin>0,.data$ymin,0))
    if(modType=="All") { #for plotting all models in 1 figure
      if(VarCalc!="None") {
      g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI, fill=.data$Source))+
        geom_line(aes(col=.data$Source))+ geom_ribbon(alpha=0.3)+xlab("Year")+
        #        geom_line(aes(y=Total.mean,col=Source),lty=2,lwd=2)+
        ylab(ytitle)}
      if(VarCalc=="None") {
      g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,fill=.data$Source))+
        geom_line(aes(col=.data$Source))+xlab("Year")+
        ylab(ytitle)
      }

    } else { #for plotting separately for each model
      if(all(is.na(yearpred$Total.mean))) {
        if(all(is.na(yearpred$Total.cv)))
          g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total))+
            geom_line()+ ylab(ytitle)  else
              g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI))+
                geom_line()+ geom_ribbon(alpha=0.3)+xlab("Year")+
                ylab(ytitle)
      } else
        g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI))+
          geom_line(aes(y=.data$Total.mean),lty=2)+
          geom_line()+ geom_ribbon(alpha=0.3)+xlab("Year")+
          ylab(ytitle)
    }
    suppressWarnings(print(g))
    if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
  }
}

#' Function to plot abundance index plus minus standard error
#'
#' @param yearpred Value
#' @param modType Value
#' @param fileName Value
#' @param subtext Value
#' @param indexVarNames Value
#' @param allVarNames Value
#' @param startYear Value
#' @param common Value
#' @param run Value
#' @param catchType Value
#' @param catchUnit Value
#' @import dplyr
#' @keywords internal
plotIndex<-function(yearpred, modType, fileName, subtext="", indexVarNames, allVarNames, startYear, common, run, catchType, catchUnit) {
  if(is.numeric(yearpred$Year) & "Year" %in% allVarNames) yearpred$Year=yearpred$Year+startYear
  if(!is.null(yearpred)) {
    if(modType %in% c("Binomial","TMBbinomial")) ytitle=paste0(common[run]," ","Positive trip index") else
      ytitle=paste0("Index ", common[run]," ",catchType[run]," (",catchUnit[run],")")
#    if(modType %in% c("Delta-Lognormal","Delta-Gamma","TMBdelta-Gamma","TMBdelta-Lognormal")) modType=paste("Delta",modType)
    yearpred<-yearpred %>% mutate(Year=as.numeric(as.character(.data$Year)),ymin=.data$Index-.data$SE,ymax=.data$Index+.data$SE) %>%
      mutate(ymin=ifelse(.data$ymin>0,.data$ymin,0))
    if(modType=="All") {
      g<-ggplot(dplyr::filter(yearpred,!.data$Source %in% c("Binomial","TMBbinomial")),aes(x=.data$Year,y=.data$Index,ymin=.data$ymin,ymax=.data$ymax,fill=.data$Source))+
        geom_line(aes(col=.data$Source))+ geom_ribbon(alpha=0.3)+xlab("Year")+
        ylab(ytitle)
    } else {
      g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Index,ymin=.data$ymin,ymax=.data$ymax))+
        geom_line()+ geom_ribbon(alpha=0.3)+xlab("Year")+
        ylab(ytitle)
    }
    if(length(indexVarNames)>1) {
      varplot=as.formula(paste0("~",paste(grep("Year",indexVarNames,invert=TRUE,value=TRUE),sep="+")))
      g=g+facet_wrap(varplot)
    }
    suppressWarnings(print(g))
    if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
  }
}

#' Function to plot total catch by all models plus a validation number
#'
#' @param yearpred Value
#' @param trueval Value
#' @param fileName Value
#' @param colName Value
#' @param allVarNames Value
#' @param startYear Value
#' @param common Value
#' @param run Value
#' @param catchType Value
#' @param catchUnit Value
#' @import dplyr
#' @keywords internal
plotSumsValidate<-function(yearpred,trueval,fileName,colName, allVarNames, startYear, common, run, catchType, catchUnit,VarCalc) {
  if(is.numeric(yearpred$Year) & "Year" %in% allVarNames)
    #yearpred$Year[yearpred$Source!="Ratio"]=yearpred$Year[yearpred$Source!="Ratio"]+startYear
    yearpred$Year[yearpred$Year<startYear]=yearpred$Year[yearpred$Year<startYear]+startYear

  yearpred<-yearpred %>%
    mutate(Year=as.numeric(as.character(.data$Year)),ymin=.data$Total-.data$Total.se,ymax=.data$Total+.data$Total.se) %>%
    mutate(ymin=ifelse(.data$ymin>0,.data$ymin,0))

  trueval<-trueval %>% rename(Total=!!colName) %>%
    mutate(ymin=NA,ymax=NA,Source="Validation",
           Total.mean=NA,TotalLCI=NA,TotalUCI=NA)
  yearpred<-bind_rows(yearpred[,c("Year","Total","Total.mean","TotalLCI","TotalUCI","ymin","ymax","Source")],
                      trueval[,c("Year","Total","Total.mean","TotalLCI","TotalUCI","ymin","ymax","Source")])

  if(VarCalc == "None"){
  g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,fill=.data$Source))+
    geom_line(aes(color=.data$Source))+
    xlab("Year")+
    ylab(paste0(common[run]," ",catchType[run]," (",catchUnit[run],")"))+
    geom_point(data=yearpred[yearpred$Source=="Validation",],aes(x=.data$Year,y=.data$Total,color=.data$Source),size=2)
  }

  if(VarCalc!="None"){
  if(all(is.na(yearpred$Total.mean))){ #if varCalc is DeltaMethod (?)
    g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI,fill=.data$Source))+
    geom_line(aes(color=.data$Source))+ geom_ribbon(alpha=0.3)+
    xlab("Year")+
    ylab(paste0(common[run]," ",catchType[run]," (",catchUnit[run],")"))+
    geom_point(data=yearpred[yearpred$Source=="Validation",],aes(x=.data$Year,y=.data$Total,color=.data$Source),size=2)
    }else{ #if VarCalc is Simulate
    g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI,fill=.data$Source))+
    geom_line(aes(color=.data$Source))+ geom_ribbon(alpha=0.3)+
    geom_line(aes(y=.data$Total.mean,color=.data$Source),lty=2)+
    xlab("Year")+
    ylab(paste0(common[run]," ",catchType[run]," (",catchUnit[run],")"))+
    geom_point(data=yearpred[yearpred$Source=="Validation",],aes(x=.data$Year,y=.data$Total,color=.data$Source),size=2)}}

  suppressWarnings(print(g))
  if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
}

#' Function to plot boxplots of RMSE and ME across folds
#'
#' @param rmse Value
#' @param me Value
#' @param fileName Value
#' @import dplyr
#' @importFrom tidyr pivot_longer
#' @keywords internal
plotCrossVal<-function(rmse,me,fileName) {
  Model <- value <- NULL
  rmse<-data.frame(rmse)%>% select_if(~!all(is.na(.)))
  me<-data.frame(me)%>% select_if(~!all(is.na(.)))
  RMSE<-tidyr::pivot_longer(rmse,everything(),names_to="Model")
  ME<-tidyr::pivot_longer(me,everything(),names_to="Model")
  df<-bind_rows(list(RMSE=RMSE,ME=ME),.id="Metric")
  g<-ggplot(df)+geom_boxplot(aes(x=Model,y=value),fill="lightgrey")+
    facet_wrap(Metric~.,ncol=1,scales="free")+
    xlab("Model")+ylab("Cross validation metrics")
  suppressWarnings(print(g))
  if(!is.null(fileName)) ggsave(fileName,height=5,width=7)
}

#' Function to fit a specified model formula and print outputs
#'
#' @param formula1 Value
#' @param formula2 Value
#' @param modType Value
#' @param obsdatval Value
#' @param outputDir Value
#' @importFrom stats update
#' @import utils
#' @keywords internal
FitModelFunc<-function(formula1,formula2,modType,obsdatval,outputDir) {
  modfit2=NULL
  formula3=update(formula2,~.+offset(log(Effort)))
  if(modType %in% c("Binomial","Delta-Lognormal","Delta-Gamma") )  {
    obsdatval$y=obsdatval$pres
    modfit1<-try(glm(formula1,data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=1000)))
  }
  if(modType %in% c("TMBBinomial","TMBdelta-Lognormal","TMBdelta-Gamma") )  {
    obsdatval$y=obsdatval$pres
    modfit1<-try(glmmTMB(formula1,data=obsdatval,family="binomial"))
  }
  if(grepl("delta",modType,ignore.case = TRUE)) {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
  }
  if(modType=="Delta-Lognormal") {
    modfit2=try(lm(formula2,data=obsdatval))
  }
  if(modType=="Delta-Gamma") {
    modfit2=try(glm(formula2,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="TMBdelta-Lognormal") {
    modfit2=try(glmmTMB(formula2,data=obsdatval))
  }
  if(modType=="TMBdelta-Gamma") {
    modfit2=try(glmmTMB(formula2,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="NegBin") {
    obsdatval$y=round(obsdatval$Catch)
    modfit1=try(glm.nb(formula3,data=obsdatval,control=glm.control(epsilon=1E-6,maxit=45),na.action=na.fail))
  }
  if(modType=="Tweedie") {
    obsdatval$y=obsdatval$cpue
    modfit1=try(cplm::cpglm(formula2,data=obsdatval))
  }
  if(modType %in% c("TMBnbinom1","TMBnbinom2") ){
    obsdatval$y=round(obsdatval$Catch)
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula3,family=TMBfamily,data=obsdatval))
  }
  if(modType =="TMBtweedie"){
    obsdatval$y=obsdatval$cpue
    TMBfamily=gsub("TMB","",modType)
    modfit1=try(glmmTMB(formula2,family=TMBfamily,data=obsdatval))
  }
  if(modType =="TMBnormal"){
    obsdatval$y=obsdatval$cpue
    modfit1=try(glmmTMB(formula2,data=obsdatval))
  }
  if(modType =="TMBlognormal"){
    obsdatval$y=log(obsdatval$cpue+0.1)
    modfit1=try(glmmTMB(formula2,data=obsdatval))
  }
  if(modType =="TMBgamma"){
    obsdatval$y=log(obsdatval$cpue+0.1)
    modfit1=try(glmmTMB(formula2,family=Gamma(link="log"),data=obsdatval))
  }
  if(class(modfit1)[1]=="try-error") modfit1=NULL
  if(class(modfit2)[1]=="try-error") modfit2=NULL
  # if(!is.null(modfit1)) {
  #   if(modType %in% c("Binomial","Delta-Lognormal","Delta-Gamma"))  #for delta models write binomial anova
  #     write.csv(anova(modfit1,test="Chi"),file=paste0(outputDir,"/BinomialAnova.csv"))
  #   if(modType %in% c("NegBin")) anova1=anova(modfit1,test="Chi")
  #   if(modType %in% c("Tweedie","TMBtweedie","TMBnbinom1","TMBnbinom2")) anova1=NULL
  #   if(modType %in% c("Delta-Lognormal","Delta-Gamma")) anova1=anova(modfit2,test="F")
  #   if(!is.null(anova1)) {
  #     write.csv(anova1,paste0(outputDir,"/anova",modType,".csv"))
  #   }
  # }
  list(modfit1=modfit1,modfit2=modfit2)
}


#' Function to make data summarizes including ratio estimate at
#'
#'stratification defined by strataVars. No pooling for missing strata
#'
#' @param obsdatval Value
#' @param logdatval Value
#' @param strataVars Value
#' @param EstimateBycatch Value
#' @param startYear Value
#' @keywords internal
MakeSummary<-function(obsdatval,logdatval,strataVars, EstimateBycatch, startYear) {
  x<-obsdatval %>%
    group_by_at(all_of(strataVars))  %>%
    summarize(OCat=sum(.data$Catch,na.rm=TRUE),
              OEff=sum(.data$Effort,na.rm=TRUE),
              OUnit=length(.data$Year),
              CPUE=mean(.data$cpue,na.rm=TRUE),
              CPse=standard.error(.data$cpue),
              Out=outlierCountFunc(.data$cpue),
              Pos=sum(.data$pres,na.rm=TRUE),
              OCatS=sd(.data$Catch,na.rm=TRUE),
              OEffS=sd(.data$Effort,na.rm=TRUE),
              Cov=cov(.data$Catch,.data$Effort)) %>%
    mutate(PFrac=.data$Pos/.data$OUnit)
  if(EstimateBycatch) {
    returnval<-logdatval  %>%
      group_by_at(all_of(strataVars)) %>%
      summarize(Eff=sum(.data$Effort,na.rm=TRUE),Units=sum(.data$SampleUnits))
    returnval<-left_join(returnval,x,by=strataVars)  %>%
      mutate(OEff=ifelse(is.na(.data$OEff),0,.data$OEff),OUnit=ifelse(is.na(.data$OUnit),0,.data$OUnit))%>%
      mutate(EFrac=.data$OEff/.data$Eff, UFrac=.data$OUnit/.data$Units) %>%
      mutate(Cat=(.data$OCat/.data$OEff)*.data$Eff,
             Cse=sqrt(ratioVar(.data$OEff,.data$Eff,.data$OUnit,.data$Units,.data$OCat/.data$OEff,.data$OEffS^2,.data$OCatS^2,.data$Cov))) %>%
      ungroup() %>% mutate(Year=as.numeric(as.character(.data$Year))) %>%
      mutate(Year=ifelse(.data$Year<startYear,.data$Year+startYear,.data$Year))
  } else {
    returnval<-x
  }
  returnval
}

#' Function to setup pooling if requested for design estimators
#'
#' @param obsdatval Value
#' @param logdatval Value
#' @param minStrataUnit Value
#' @param designVars Value
#' @param pooledVar Value
#' @param poolTypes Value
#' @param adjacentNum Value
#' @keywords internal
getPooling<-function(obsdatval,logdatval,minStrataUnit,designVars,
  pooledVar,poolTypes,adjacentNum) {
 obsdatval<-data.frame(obsdatval)
 logdatval<-data.frame(logdatval)
 poolingVars<-c(designVars,pooledVar[!is.na(pooledVar) &!pooledVar %in% designVars])
 if("Year" %in% poolingVars) {
   if(is.factor(obsdatval$Year)) yearFactor<-TRUE else yearFactor<-FALSE
   if(is.factor(obsdatval$Year)) obsdatval$Year=as.numeric(as.character(obsdatval$Year))
   if(is.factor(logdatval$Year)) logdatval$Year=as.numeric(as.character(logdatval$Year))
 }
 poolingSum<-logdatval %>% group_by_at(all_of(poolingVars)) %>%
  summarize(totalUnits=sum(.data$SampleUnits),totalEffort=sum(.data$Effort))
 x<-obsdatval %>%group_by_at(all_of(poolingVars)) %>%
   summarize(units=n(),effort=sum(.data$Effort))
 poolingSum<-left_join(poolingSum,x,by=poolingVars) %>%
   mutate(units=replace_na(.data$units,0),effort=replace_na(.data$effort,0)) %>%
   mutate(needs.pooling=ifelse(.data$units>minStrataUnit, FALSE,TRUE),
          pooled.n=ifelse(.data$units>minStrataUnit, units,NA),
          poolnum=NA,pooledTotalUnits=NA,pooledTotalEffort=NA) %>%
   ungroup()
  poolingSum<-as.data.frame(poolingSum)
  poolingSum$poolnum[!poolingSum$needs.pooling]<-0
  includePool<-list()
  for(i in which(!poolingSum$needs.pooling))  {
   poolingSum$pooledTotalUnits[i]<-poolingSum$totalUnits[i]
   poolingSum$pooledTotalEffort[i]<-poolingSum$totalEffort[i]
   bb<-1:nrow(obsdatval)
   for(j in 1:length(designVars)) {
     bb<-bb[bb %in% which(obsdatval[,designVars[j]]==poolingSum[i,designVars[j]])]
   }
   includePool[[i]]<-obsdatval[bb,]
  }
  for(vari in 1:length(designVars))  {
   keepVars<-designVars[(1:length(designVars))>vari]
   for(i in which(poolingSum$needs.pooling))  {
     if(poolTypes[1]=="none") {
       aa<-which(poolingSum[,designVars[1]] == poolingSum[i,designVars[1]])
       bb<-which(obsdatval[,designVars[1]] == poolingSum[i,designVars[1]])
     }
     if(poolTypes[1]=="all") {
      aa<-1:nrow(poolingSum)
      bb<-1:nrow(obsdatval)
    }
    if(poolTypes[1]=="pooledVar") {
      aa<-which(poolingSum[,pooledVar[1]]==poolingSum[i,pooledVar[1]])
      bb<-which(obsdatval[,pooledVar[1]]==poolingSum[i,pooledVar[1]])
    }
    if(poolTypes[1]=="adjacent") {
      aa<-which(poolingSum[,designVars[1]] >= poolingSum[i,designVars[1]]-adjacentNum[1] &
                  poolingSum[,designVars[1]] <= poolingSum[i,designVars[1]]+adjacentNum[1])
      bb<-which(obsdatval[,designVars[1]] >= poolingSum[i,designVars[1]]-adjacentNum[1] &
                                  obsdatval[,designVars[1]] <= poolingSum[i,designVars[1]]+adjacentNum[1])
    }
    if(vari>1) {
      for(var2 in 2:vari) {
        if(poolTypes[var2]=="none") {
          aa<-aa[aa %in% which(poolingSum[,designVars[var2]] == poolingSum[i,designVars[var2]])]
          bb<-bb[bb %in% which(obsdatval[,designVars[var2]] == poolingSum[i,designVars[var2]])]
        }
        if(poolTypes[var2]=="pooledVar") {
          aa<-aa[aa %in% which(poolingSum[,pooledVar[var2]]==poolingSum[i,pooledVar[var2]])]
          bb<-bb[bb %in% which(obsdatval[,pooledVar[var2]]==poolingSum[i,pooledVar[var2]])]
        }
        if(poolTypes[var2]=="adjacent") {
          aa<-aa[aa %in% which(poolingSum[,designVars[var2]] >= poolingSum[i,designVars[var2]]-adjacentNum[var2] &
                      poolingSum[,designVars[var2]] <= poolingSum[i,designVars[var2]]+adjacentNum[var2])]
          bb<-bb[bb %in% which(obsdatval[,designVars[var2]] >= poolingSum[i,designVars[var2]]-adjacentNum[var2] &
                      obsdatval[,designVars[var2]] <= poolingSum[i,designVars[var2]]+adjacentNum[var2])]
        }
      }
    }
    if(length(keepVars)>0) {
    for(j in 1:length(keepVars)) {
      aa<-aa[aa %in% which(poolingSum[,keepVars[j]]==poolingSum[i,keepVars[j]])]
      bb<-bb[bb %in% which(obsdatval[,keepVars[j]]==poolingSum[i,keepVars[j]])]
    }}
    if(length(bb)>0) {
     includePool[[i]]<-obsdatval[bb,]
     poolingSum$pooled.n[i]<-nrow(includePool[[i]])
     poolingSum$needs.pooling[i]<-ifelse(poolingSum$pooled.n[i]>=minStrataUnit,FALSE,TRUE)
     poolingSum$pooledTotalEffort[i]<-sum(poolingSum$totalEffort[aa])
     poolingSum$pooledTotalUnits[i]<-sum(poolingSum$totalUnits[aa])
    } else poolingSum$needs.pooling[i]<-TRUE
   }
   poolingSum$poolnum[!poolingSum$needs.pooling &is.na(poolingSum$poolnum)]<-vari
  }
  includePool<-bind_rows(includePool,.id="stratum")
  poolingSum$stratum<-1:nrow(poolingSum)
  if("Year" %in% poolingVars) {
    if(yearFactor) {
      poolingSum$Year<-factor(poolingSum$Year)
      includePool$Year<-factor(includePool$Year)
    }
  }
  list(poolingSum,includePool)
}


#' Function to make design based estimates of bycatch from the
#' ratio estimator of Pennington Delta estimator, pooling as
#' needed for strata missing data.
#' stratification defined by designVars, then aggregated to strataVars
#'
#' @param obsdatval Value
#' @param logdatval Value
#' @param strataVars Value
#' @param designVars Value
#' @param designPooling Value
#' @param minStrataUnit Value
#' @param startYear Value
#' @param poolingSum Value
#' @param includePool Value
#' @keywords internal
getDesignEstimates<-function(obsdatval,logdatval,strataVars,designVars=NULL,
                             designPooling,minStrataUnit=1,startYear,
                             poolingSum=NULL,includePool=NULL) {
  if(!designPooling) {
    x<-obsdatval %>%
      group_by_at(all_of(designVars))  %>%
      summarize(OCat=sum(.data$Catch,na.rm=TRUE),
                OEff=sum(.data$Effort,na.rm=TRUE),
                OUnit=length(.data$Year),
                CPUE=mean(.data$cpue,na.rm=TRUE),
                CPse=standard.error(.data$cpue),
                Pos=sum(.data$pres,na.rm=TRUE),
                OCatS=sd(.data$Catch,na.rm=TRUE),
                OEffS=sd(.data$Effort,na.rm=TRUE),
                Cov=cov(.data$Catch,.data$Effort, use="complete.obs" ),
                deltaMeanCPUE=deltaEstimatorMean(.data$cpue),
                deltaVar=deltaEstimatorVar(.data$cpue),
                deltaSE2=deltaEstimatorSE2(.data$cpue)) %>%
      mutate(PFrac=.data$Pos/.data$OUnit)
    returnval<-logdatval  %>%
      group_by_at(all_of(designVars)) %>%
      summarize(Eff=sum(.data$Effort,na.rm=TRUE),
                Units=sum(.data$SampleUnits))
    returnval<-left_join(returnval,x,by=designVars)  %>%
      mutate(OEff=ifelse(is.na(.data$OEff),0,.data$OEff),
             OUnit=ifelse(is.na(.data$OUnit),0,.data$OUnit)) %>%
      mutate(ratioMean=(.data$OCat/.data$OEff)*.data$Eff,
             ratioSE=sqrt(ratioVar(.data$OEff,.data$Eff,.data$OUnit,.data$Units,
                                   .data$OCat/.data$OEff,.data$OEffS^2,.data$OCatS^2,.data$Cov)),
             deltaMean=.data$deltaMeanCPUE*.data$Eff,
             deltaSE=sqrt(.data$deltaSE2)*.data$Eff)
    returnval<-replace_na(returnval,list(ratioMean=0,ratioSE=0,deltaMean=0,deltaSE=0))
  } else {  #For pooling
  poolVars<-designVars
  logdatval<-left_join(logdatval,poolingSum[,c(designVars,"stratum")],by=designVars)
  x<-obsdatval %>%
      group_by_at(all_of(poolVars) ) %>%
      summarize(OCat=sum(.data$Catch,na.rm=TRUE),
                OEff=sum(.data$Effort,na.rm=TRUE),
                OUnit=length(.data$Year),
                CPUE=mean(.data$cpue,na.rm=TRUE),
                CPse=standard.error(.data$cpue),
                Pos=sum(.data$pres,na.rm=TRUE),
                OCatS=sd(.data$Catch,na.rm=TRUE),
                OEffS=sd(.data$Effort,na.rm=TRUE),
                Cov=cov(.data$Catch,.data$Effort, use="complete.obs" )) %>%
      mutate(PFrac=.data$Pos/.data$OUnit)
     y<-includePool %>%
      group_by(.data$stratum) %>%
      summarize(deltaMeanCPUE=deltaEstimatorMean(.data$cpue),
                deltaVar=deltaEstimatorVar(.data$cpue),
                deltaSE2=deltaEstimatorSE2(.data$cpue),
                pCat=sum(.data$Catch,na.rm=TRUE),
                pEff=sum(.data$Effort,na.rm=TRUE),
                pUnit=n(), #eab
                pCatS=sd(.data$Catch,na.rm=TRUE),
                pEffS=sd(.data$Effort,na.rm=TRUE),
                pCov=cov(.data$Catch,.data$Effort, use="complete.obs" )) %>%
      ungroup() %>%
      mutate(stratum=as.numeric(as.character(.data$stratum)))
    poolVals<-logdatval  %>%
      group_by_at(all_of(c(poolVars,"stratum"))) %>%
      summarize(Eff=sum(.data$Effort,na.rm=TRUE),
                Units=sum(.data$SampleUnits))
    poolVals<-left_join(poolVals,x,by=poolVars)
    poolVals<-left_join(poolVals,poolingSum[,c(poolVars,"pooledTotalEffort","pooledTotalUnits")],by=poolVars)
    poolVals<-left_join(poolVals,y,by="stratum")
    poolVals<-poolVals %>%
      mutate(OEff=ifelse(is.na(.data$OEff),0,.data$OEff),
             OUnit=ifelse(is.na(.data$OUnit),0,.data$OUnit)) %>%
      mutate(ratioMean=(.data$pCat/.data$pEff)*.data$Eff,
             ratioSE=sqrt(ratioVar(.data$pEff,.data$pooledTotalEffort,.data$pUnit,.data$pooledTotalUnits,
                                   .data$pCat/.data$pEff,.data$pEffS^2,.data$pCatS^2,.data$pCov))*
                             .data$Eff/.data$pooledTotalEffort,
             deltaMean=.data$deltaMeanCPUE*.data$Eff,
             deltaSE=sqrt(.data$deltaSE2)*.data$Eff)
    returnval<-replace_na(poolVals,list(ratioMean=0,ratioSE=0,deltaMean=0,deltaSE=0))
  }
  if(all(is.na(strataVars))) {
    returnval<-returnval %>%
      ungroup() %>%
      summarize(ratioMean=sum(.data$ratioMean,na.rm=TRUE),
                ratioSE=sqrt(sum(.data$ratioSE^2,na.rm=TRUE)),
                deltaMean=sum(.data$deltaMean,na.rm=TRUE),
                deltaSE=sqrt(sum(.data$deltaSE^2,na.rm=TRUE)))
  } else {
  returnval<-returnval %>%
    group_by_at(all_of(strataVars)) %>%
    summarize(ratioMean=sum(.data$ratioMean,na.rm=TRUE),
              ratioSE=sqrt(sum(.data$ratioSE^2,na.rm=TRUE)),
              deltaMean=sum(.data$deltaMean,na.rm=TRUE),
              deltaSE=sqrt(sum(.data$deltaSE^2,na.rm=TRUE))) %>%
    ungroup()
  }
    if("Year" %in% strataVars) {
      returnval<-returnval %>%
        mutate(Year=as.numeric(as.character(.data$Year))) %>%
        mutate(Year=ifelse(.data$Year<startYear,.data$Year+startYear,.data$Year))
  }
  returnval
}




#' Calculate variance of ratio estimator, from data already summarized by strata
#'
#' Variables are x=sum(obs Effort),X=sum(log Effort), n=observed sample units. N=log sample units, Rhat=mean(obs Catch)/mean(obs Effort), sx2, sy2, and sxy are the observed variance in effort and catch, and covariance
#'
#' @param x Value
#' @param X Value
#' @param n Value
#' @param N Value
#' @param Rhat Value
#' @param sx2 Value
#' @param sy2 Value
#' @param sxy Value
#' @keywords internal
ratioVar<-function(x,X,n,N,Rhat,sx2,sy2,sxy) {
  X^2*(1-n/N)/(x^2/n)*(sy2+Rhat^2*sx2-2*Rhat*sxy)
}

#' Return a column of R squared values from input of a MuMin dredge table
#'
#' @param dredgeTable dredge output from MuMin
#' @param obsdatval observer data
#' @param funcName function used in fitting (e.g. glm or glmmTMB)
#' @keywords internal
addR2<-function(dredgeTable,obsdatval,funcName) {
  R2<-NA
  for(i in 1:nrow(dredgeTable)) {
    if(funcName=="cpglm") {
      R2[i]=NA
    } else {
      mod1<-get.models(dredgeTable,subset=i)[[1]]
      mod2<-try(do.call(funcName,args=list(formula= formula(mod1$call),
                                           data=obsdatval)))
      if(funcName=="glmmTMB") {
        if(ncol(model.matrix(mod2))==length(fixef(mod2)[[1]]))
          R2[i]=r.squaredGLMM(mod2)[1,"R2c"] else
            R2[i]=NA

      } else
        R2[i]=r.squaredGLMM(mod2)[1,"R2c"]
    }
  }
  R2
}


#' Return a table of model parameters and summaries
#'
#' @param modfits A list of fitted model objects
#' @param modTypes A corresponding vector of the model types as specified in modelTry
#' @keywords internal
getModelSummaryTable<-function(modfits,modTypes) {
  modSum<-data.frame(Model=modTypes,scale=NA,phi=NA,loglike=NA,df.resid=NA)
  for(i in 1:length(modTypes)) {
    mod1<-modfits[[i]]
    modType<-modTypes[i]
    modSum$loglike[i]<-as.vector(logLik(mod1))
    modSum$df.resid[i]<-df.residual(mod1)
    if(!grepl("binomial",modType,ignore.case = TRUE))
      modSum$scale[i]<-sigma(mod1)
    if(modType == "Tweedie") {
      modSum$scale[i]=mod1$p
      modSum$phi[i]=mod1$phi
    }
    if(modType == "NegBin" )  {
      modSum$scale[i]<-mod1$theta
    }
    if(modType == "TMBtweedie" )  {
      modSum$phi[i]<-glmmTMB::family_params(mod1)
    }
  }
  modSum
}

#' Read in all the R objects created during runs of the bycatchEstimator for further analysis
#'
#' @param baseDir The base directory for the runs, same as bycatchSetup.
#' @param runName The run name, same as bycatchSetup.
#' @param runDate  The date when the model was run. Defualts to current date, but can be set to read in models previously run.
#' @param loadDesign  TRUE to read in design-based estimator results.
#' @param loadModel TRUE to read in model-based estimator results.
#' @keywords internal
loadOutputs<-function(baseDir = getwd(),
                      runName,
                      runDate =  Sys.Date(),
                      loadDesign = TRUE,
                      designScenario = NULL,
                      loadModel = TRUE,
                      modelScenario = NULL) {
  #Check that setup file exists and read in.
  outDir<-paste0(baseDir, paste("/Output", runName))
  if(!dir.exists(outDir)) stop(paste("Directory",outDir,"not found."))
  setupFile<-paste0(runDate,"_BycatchSetupSpecification.rds")
  if(!file.exists(paste0(outDir,"/",setupFile))) stop(paste("Setup file",setupFile ,"not found in",outDir,"."))
  setupObj<-readRDS(file=paste0(outDir,"/",Sys.Date(),"_BycatchSetupSpecification.rds"))
  list2env(setupObj$bycatchInputs, envir = .GlobalEnv)
  list2env(setupObj$bycatchOutputs, envir = .GlobalEnv)
  #If doing design based, check that design file exists and read in.
  if(loadDesign) {
    designFile<-paste0(runDate,"_BycatchDesign",designScnenario,".rds")
    if(!file.exists(paste0(outDir,"/",designFile))) stop(paste("Design file",designFile ,"not found in",outDir,"."))
    designObj<-readRDS(file=paste0(outDir,"/",Sys.Date(),"_BycatchDesign",designFile,".rds"))
    list2env(designObj$designInputs, envir = .GlobalEnv)
    list2env(designObj$designOutputs, envir = .GlobalEnv)
  }
  if(loadModel) {
    modelFile<-paste0(runDate,"_BycatchFit",modelScenario,".rds")
    if(!file.exists(paste0(outDir,"/",modelFile))) stop(paste("Model file",modelFile ,"not found in",outDir,"."))
    modelObj<-readRDS(file=paste0(outDir,"/",Sys.Date(),"_BycatchFit",modelScenario,".rds"))
    list2env(modelObj$modelInputs, envir = .GlobalEnv)
    list2env(modelObj$modelOutputs, envir = .GlobalEnv)
  }
}
