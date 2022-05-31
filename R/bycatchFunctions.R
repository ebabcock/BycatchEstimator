
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
#' @param useParallel Value
#' @param selectCriteria Value
#' @param varExclude Value
#' @param printOutput Value
#' @param catchType Value
#' @param common Value
#' @param dirname Value
#' @param run Value
#' @import MuMIn parallel doParallel tweedie glmmTMB
#' @importFrom reshape2 colsplit
#' @importFrom stats anova na.fail as.formula coef Gamma glm.control formula lm glm vcov
#' @importFrom MASS glm.nb
#' @keywords internal
findBestModelFunc<-function(obsdatval, modType, requiredVarNames, allVarNames, complexModel, useParallel, selectCriteria, varExclude, printOutput=FALSE, catchType = NULL, common = NULL, dirname = NULL, run = NULL) {

  offset<-TMBfamily<-NULL
  keepVars=requiredVarNames
  extras=c("AICc","AIC", "BIC")

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

  if(modType %in% c("NegBin","TMBnbinom1","TMBnbinom2"))
    obsdatval$y=round(obsdatval$Catch)
  if(modType %in% c("Tweedie","TMBtweedie","Normal"))
    obsdatval$y=obsdatval$cpue
  if(modType =="Binomial")   {
    obsdatval$y=obsdatval$pres
    funcName="glm"
    args=list(formula="",data=obsdatval,family="binomial",control=list(epsilon = 1e-6,maxit=100),na.action=na.fail)
  }
  if(modType=="Normal") {
    funcName="lm"
    args=list(formula="",data=obsdatval,na.action=na.fail)
  }
  if(modType=="Lognormal") {
    obsdatval$y=log(obsdatval$cpue+0.1)
    funcName="lm"
    args=list(formula="",data=obsdatval,na.action=na.fail)
  }
  if(modType=="Gamma") {
    obsdatval$y=obsdatval$cpue+0.1
    funcName="glm"
    args=list(formula="",data=obsdatval,family=Gamma(link="log"),na.action=na.fail)
  }
  if(modType=="Delta-Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    funcName="lm"
    args=list(formula="",data=obsdatval,na.action=na.fail)
  }
  if(modType=="Delta-Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    funcName="glm"
    args=list(formula="",data=obsdatval,family=Gamma(link="log"),na.action=na.fail)
  }
  if(modType == "NegBin") {
    funcName="glm.nb"
    offset="+offset(log(Effort))"
    keepVars=c(requiredVarNames,"offset(log(Effort))")
    args=list(formula="",data=obsdatval,control=glm.control(epsilon=1E-6,maxit=30),na.action=na.fail)
  }
  if(modType=="Tweedie") {
    funcName="cpglm"
    args=list(formula="",data=obsdatval,na.action=na.fail)
  }
  if(modType %in% c("TMBnbinom1","TMBnbinom2") ){
    funcName="glmmTMB"
    TMBfamily=gsub("TMB","",modType)
    offset="+offset(log(Effort))"
    keepVars=paste0("cond(",c(requiredVarNames,"offset(log(Effort))"),")")
    args=list(formula="",family=TMBfamily,data=obsdatval,na.action=na.fail)
  }
  if(modType %in% c("TMBtweedie") ){
    funcName="glmmTMB"
    TMBfamily=gsub("TMB","",modType)
    keepVars=paste0("cond(",requiredVarNames,")")
    args=list(formula="",family=TMBfamily,data=obsdatval,na.action=na.fail)
  }
  formulaList<-list(as.formula(paste("y~",paste(getAllTerms(complexModel),collapse="+"),offset)),
                    as.formula(paste("y~",paste(allVarNames,collapse="+"),offset)),
                    as.formula(paste("y~",paste(allVarNames[!allVarNames %in% varExclude],collapse="+"),offset)),
                    as.formula(paste("y~",paste(requiredVarNames,collapse="+"),offset)), NA)
  args$formula=formulaList[[1]]
  modfit1<-try(do.call(funcName,args))
  for(i in 2:(length(formulaList))-1) {
    if(class(modfit1)[1] %in% c("glm","lm","glm.nb")) {
      if(modfit1$rank<length(coef(modfit1))) class(modfit1)<-"try-error"
    }
    if(class(modfit1)[1] == "cpglm")  {
      if(length(coef(modfit1))!=dim(modfit1$vcov)[1])  class(modfit1)<-"try-error"
    }
    if(class(modfit1)[1]=="try-error")   {
      args$formula=formulaList[[i]]
      modfit1<-try(do.call(funcName,args))
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
    if(modType=="Tweedie") modfit1<-cpglm(formula(modfit1),data=obsdatval,na.action=na.fail)
    if(modType %in% c("TMBnbinom1","TMBnbinom2","TMBtweedie") )
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
    if(printOutput & !is.null(modfit2)) {
      write.csv(modfit2,paste0(dirname[[run]],common[run],catchType[run],"ModelSelection",modType,".csv"))
      if(modType %in% c("Binomial","NegBin")) anova1=anova(modfit3,test="Chi")
      if(modType %in% c("Tweedie","TMBnbinom1","TMBnbinom2","TMBtweedie")) anova1=NULL
      if(modType %in% c("Normal","Lognormal","Gamma","Delta-Lognormal","Delta-Gamma")) anova1=anova(modfit3,test="F")
      if(!is.null(anova1)) {
        write.csv(anova1,paste0(dirname[[run]],common[run],catchType[run],modType,"Anova.csv"))
      }
    }
    returnval=list(modfit3,modfit2)
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
#' @import tidyr
#' @importFrom stats predict model.matrix rbinom sigma rnorm rlnorm rnbinom quantile
#' @importFrom MASS mvrnorm gamma.shape
#' @keywords internal
makePredictionsSimVarBig<-function(modfit1, modfit2=NULL, newdat, modtype, obsdatval, includeObsCatch, nsim, requiredVarNames, CIval, printOutput=TRUE, catchType, common, dirname, run) {
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
  if(modtype=="Tweedie" ) response1<-data.frame(cplm::predict(modfit1,newdata=newdat,type="response",se.fit=TRUE)) else
   response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(dim(response1)[2]==1) {
      names(response1)="fit"
      if(modtype=="Tweedie")
        response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp",offsetval=NULL, nsim=nsim) else
          response1$se.fit=rep(NA,dim(response1)[2])
    }
    if(!is.null(modfit2))  {
      response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
      names(response2)=paste0(names(response2),"2")
    }
 if(!any(is.na(response1$se.fit)) & !max(response1$se.fit/response1$fit,na.rm=TRUE)>10000)  {
  #Set up model matrices for simulation
  yvar=sub( " ", " ",formula(modfit1) )[2]
  newdat<-cbind(y=rep(1,nObs),newdat)
  names(newdat)[1]=yvar
  a=model.matrix(formula(modfit1),data=newdat)
  if(!is.null(modfit2)) {
   yvar=sub( " ", " ",formula(modfit2) )[2]
   if(! yvar %in% names(newdat)) {
    newdat<-cbind(y=rep(1,nObs),newdat)
    names(newdat)[1]=yvar
   }
   b=model.matrix(formula(modfit2),data=newdat)
  }
  #Get predictions sim
  if(modtype == "Binomial") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=.data$fit,TotalVar=.data$se.fit^2+.data$fit*(1-.data$fit))
      sim=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )
     }
  if(modtype=="Normal") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*.data$fit,
               TotalVar=.data$Effort^2*(.data$se.fit^2+sigma(modfit1)^2))
      sim=replicate(nsim,rnorm(nObs,
        mean=as.vector(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))),
        sd=sigma(modfit1)))*newdat$Effort
  }
  if(modtype=="Lognormal") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*(lnorm.mean(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))-0.1),
               TotalVar=.data$Effort^2*lnorm.se(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))^2)
      sim=replicate(nsim,rlnorm(nObs,
                                meanlog=as.vector((a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))),
                                sdlog=sigma(modfit1))-0.1)*newdat$Effort
  }
  if(modtype=="Gamma") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*(.data$fit-0.1),
               TotalVar=.data$Effort^2*(.data$se.fit^2+.data$fit*gamma.shape(modfit1)[[1]]))
      sim=replicate(nsim,newdat$Effort*(simulateGammaDraw(modfit1,nObs,a)-0.1) )
  }
  if(modtype == "Delta-Lognormal") {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue=lnorm.mean(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)),
               pos.cpue.se=lnorm.se(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)),
               prob.se=sqrt(.data$se.fit^2+.data$fit*(1-.data$fit))) %>%
        mutate(Total=.data$Effort*.data$fit*.data$pos.cpue,
               TotalVar=.data$Effort^2*lo.se(.data$fit,.data$prob.se,.data$pos.cpue,.data$pos.cpue.se)^2)
      sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )
      sim2=replicate(nsim,newdat$Effort*exp(rnorm(nObs,b %*%
           mvrnorm(1,coef(modfit2),vcov(modfit2)),sigma(modfit2)) ) )
      sim=sim1*sim2
  }
  if(modtype == "Delta-Gamma") {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue.se=sqrt(.data$se.fit2^2+.data$fit2*gamma.shape(modfit2)[[1]]),
               prob.se=sqrt(.data$se.fit^2+.data$fit*(1-.data$fit))) %>%
        mutate(Total=.data$Effort*.data$fit*.data$fit2,
               TotalVar=.data$Effort^2*lo.se(.data$fit,.data$prob.se,.data$fit2,.data$pos.cpue.se)^2)
      sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )
      sim2=replicate(nsim,newdat$Effort*simulateGammaDraw(modfit2,nObs,b) )
      sim=sim1*sim2
  }
  if(modtype=="NegBin") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$fit,TotalVar=.data$se.fit^2+.data$fit+.data$fit^2/modfit1$theta)
      sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))*newdat$Effort,
        size=modfit1$theta))  #Simulate negative binomial data
  }
  if(modtype=="Tweedie") {
      allpred=cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*.data$fit,
               TotalVar=.data$Effort^2*(.data$se.fit^2+modfit1$phi*.data$fit^modfit1$p))
       sim=replicate(nsim,rtweedie(nObs,power=modfit1$p,
        mu=as.vector(exp(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))),
         phi=modfit1$phi))*newdat$Effort
  }
  if(modtype=="TMBnbinom1") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$fit,
               TotalVar=.data$se.fit^2+.data$fit+.data$fit*sigma(modfit1))
      sim = replicate(nsim,simulateNegBin1Draw(modfit1,nObs,a,newdat$Effort))
  }
  if(modtype=="TMBnbinom2") {
       allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$fit,
               TotalVar=.data$se.fit^2+.data$fit+.data$fit^2/sigma(modfit1))
      sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,fixef(modfit1)[[1]],
        vcov(modfit1)[[1]]))*newdat$Effort, size=sigma(modfit1)))
  }
  if(modtype=="TMBtweedie") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$Effort*.data$fit,
               TotalVar=.data$Effort^2*(.data$se.fit^2+sigma(modfit1)*.data$fit^(glmmTMB:::.tweedie_power(modfit1))))
       sim=replicate(nsim, simulateTMBTweedieDraw(modfit1,nObs,a,newdat$Effort) )
  }
    if(includeObsCatch & modtype!="Binomial") {
      obsdatvalyear=obsdatval[obsdatval$Year==years[i],]
      d=match(allpred$matchColumn,obsdatvalyear$matchColumn)
      allpred$Total[!is.na(d)]= allpred$Total[!is.na(d)] + obsdatvalyear$Catch[d[!is.na(d)]]
      sim[!is.na(d),]= sim[!is.na(d),] + obsdatvalyear$Catch[d[!is.na(d)]]
    }
    if(includeObsCatch & modtype=="Binomial") {
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
       write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modtype,"StratumSummary.csv"))
       write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modtype,"AnnualSummary.csv"))
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
#' @importFrom stats delete.response terms qnorm
#' @keywords internal
makePredictionsDeltaVar<-function(modfit1, newdat, modtype,  obsdatval, includeObsCatch, requiredVarNames, CIval, printOutput=TRUE, catchType, common, dirname, run) {
  if(modtype %in% c("Delta-Lognormal","Delta-Gamma")) stop("No delta-method variance available")
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
     cplm::predvallink = predict(modfit1,newdat=newdat)
     cplm::predval = predict(modfit1,newdat=newdat,type="response")

    } else {
     predvallink = predict(modfit1,newdat=newdat)
     predval = predict(modfit1,newdat=newdat,type="response")
    }
    if(modtype %in% c("TMBtweedie","TMBnbinom1","TMBnbinom2"))
      vcovval = a %*% vcov(modfit1)[[1]] %*% t(a) else
        vcovval = a %*% vcov(modfit1) %*% t(a)
    if(modtype == "Binomial") {
      residvar =  predval * (1-predval)  #Binomial variance
      deriv =  as.vector(exp(predvallink)/(exp(predvallink)+1)^2)
    }
    if(modtype == "Normal") {
      predval=predval*newdat$Effort
      residvar =  rep(sigma(modfit1)^2,nrow(newdat))*newdat$Effort^2
      deriv =  rep(1,nrow(newdat))
    }
    if(modtype == "Lognormal" ) {
      temp = predict(modfit1,newdata=newdat,se.fit=TRUE)
      predval = (lnorm.mean(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2))-0.1)*newdat$Effort
      # residvar = lnorm.se(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2))^2*newdat$Effort^2
      # deriv =  exp(temp$fit)*newdat$Effort
      deriv = (lnorm.mean(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2)))*newdat$Effort
      residvar = lnorm.se(temp$fit,sqrt(temp$se.fit^2+sigma(modfit1)^2))^2*newdat$Effort^2
    }
    if(modtype == "Gamma" ) {
      predval = (predval-0.1) * newdat$Effort
      predval[predval<0]<-0
      residvar =exp(predvallink)*gamma.shape(modfit1)[[1]]*newdat$Effort^2
      deriv =  exp(predvallink)*newdat$Effort
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
      residvar =  sigma(modfit1)*predval^(glmmTMB:::.tweedie_power(modfit1))*newdat$Effort^2
      predval = predval * newdat$Effort
      deriv =  predval  #derivative of exp(x) is exp(x)
    }
    if(includeObsCatch & modtype!="Binomial") {
      obsdatvalyear=obsdatval[obsdatval$Year==years[i],]
      d=match(allpred$matchColumn,obsdatvalyear$matchColumn)
      #d=match(logdat$matchColumn,obsdatvalyear$matchColumn)
      d=a[!is.na(d)]
      allpred$Total[a]= allpred$Total[d] + obsdatvalyear$Catch
    }
    if(includeObsCatch & modtype=="Binomial") {
      obsdatvalyear=obsdatval[obsdatval$Year==years[i],]
      d=match(allpred$matchColumn,obsdatvalyear$matchColumn)
      #d=match(logdat$matchColumn,obsdatvalyear$matchColumn)
      d=d[!is.na(d)]
      allpred$Total[d]= obsdatvalyear$pres
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
    if(printOutput)  write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modtype,"StratumSummary.csv"))

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
  if(printOutput) {
    write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modtype,"AnnualSummary.csv"))
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
#' @import tidyr
#' @importFrom stats predict model.matrix rbinom sigma rnorm rlnorm rnbinom quantile
#' @importFrom MASS mvrnorm gamma.shape
#' @keywords internal
makePredictionsSimVar<-function(modfit1, modfit2=NULL, modtype, newdat, obsdatval=NULL, includeObsCatch, nsim, requiredVarNames, CIval, printOutput=TRUE, catchType, common, dirname, run) {
  #Separate out sample units
  if(includeObsCatch)    newdat$Effort=newdat$unsampledEffort/newdat$SampleUnits else
    newdat$Effort=newdat$Effort/newdat$SampleUnits
  newdat=uncount(newdat,.data$SampleUnits)
  nObs=dim(newdat)[1]
  #Get predictions
  if(modtype=="Tweedie")  response1<-data.frame(cplm::predict(modfit1,newdata=newdat,type="response",se.fit=TRUE)) else
  response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(dim(response1)[2]==1) {
      names(response1)="fit"
      if(modtype=="Tweedie")
        response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp",offsetval=NULL, nsim=nsim) else
          response1$se.fit=rep(NA,dim(response1)[2])
    }
    if(!is.null(modfit2))  {
      response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
      names(response2)=paste0(names(response2),"2")
    }
if(!any(is.na(response1$se.fit)) & !max(response1$se.fit/response1$fit,na.rm=TRUE)>10000)  {
  #Set up model matrices for simulation
  yvar=sub( " ", " ",formula(modfit1) )[2]
  newdat<-cbind(y=rep(1,nObs),newdat)
  names(newdat)[1]=yvar
  a=model.matrix(formula(modfit1),data=newdat)
  if(!is.null(modfit2)) {
   yvar=sub( " ", " ",formula(modfit2) )[2]
   if(! yvar %in% names(newdat)) {
    newdat<-cbind(y=rep(1,nObs),newdat)
    names(newdat)[1]=yvar
   }
   b=model.matrix(formula(modfit2),data=newdat)
  }
  #Get predictions sim
  if(modtype == "Binomial") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=.data$fit,TotalVar=.data$se.fit^2+.data$fit*(1-.data$fit))
      sim=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )
     }
  if(modtype=="Normal") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*.data$fit,
               TotalVar=.data$Effort^2*(.data$se.fit^2+sigma(.data$modfit1)^2))
      sim=replicate(nsim,rnorm(nObs,
        mean=as.vector(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))),
        sd=sigma(modfit1)))*newdat$Effort
  }
  if(modtype=="Lognormal") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*(lnorm.mean(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))-0.1),
               TotalVar=.data$Effort^2*lnorm.se(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))^2)
      sim=replicate(nsim,rlnorm(nObs,
                                meanlog=as.vector((a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))),
                                sdlog=sigma(modfit1))-0.1)*newdat$Effort
  }
  if(modtype=="Gamma") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*(.data$fit-0.1),
               TotalVar=.data$Effort^2*(.data$se.fit^2+.data$fit*gamma.shape(modfit1)[[1]]))
      sim=replicate(nsim,newdat$Effort*(simulateGammaDraw(modfit1,nObs,a)-0.1) )
  }
  if(modtype == "Delta-Lognormal") {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue=lnorm.mean(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)),
               pos.cpue.se=lnorm.se(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)),
               prob.se=sqrt(.data$se.fit^2+.data$fit*(1-.data$fit))) %>%
        mutate(Total=.data$Effort*.data$fit*.data$pos.cpue,
               TotalVar=.data$Effort^2*lo.se(.data$fit,.data$prob.se,.data$pos.cpue,.data$pos.cpue.se)^2)
      sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )
      sim2=replicate(nsim,newdat$Effort*exp(rnorm(nObs,b %*%
           mvrnorm(1,coef(modfit2),vcov(modfit2)),sigma(modfit2)) ) )
      sim=sim1*sim2
  }
  if(modtype == "Delta-Gamma") {
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue.se=sqrt(.data$se.fit2^2+.data$fit2*gamma.shape(modfit2)[[1]]),
               prob.se=sqrt(.data$se.fit^2+.data$fit*(1-.data$fit))) %>%
        mutate(Total=.data$Effort*.data$fit*.data$fit2,
               TotalVar=.data$Effort^2*lo.se(.data$fit,.data$prob.se,.data$fit2,.data$pos.cpue.se)^2)
      sim1=replicate(nsim,rbinom(nObs,1,ilogit(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1))) ) )
      sim2=replicate(nsim,newdat$Effort*simulateGammaDraw(modfit2,nObs,b) )
      sim=sim1*sim2
  }
  if(modtype=="NegBin") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$fit,TotalVar=.data$se.fit^2+.data$fit+.data$fit^2/modfit1$theta)
      sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))*newdat$Effort,
        size=modfit1$theta))  #Simulate negative binomial data
  }
  if(modtype=="Tweedie") {
      allpred=cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*.data$fit,
               TotalVar=.data$Effort^2*(.data$se.fit^2+modfit1$phi*.data$fit^modfit1$p))
       sim=replicate(nsim,rtweedie(nObs,power=modfit1$p,
        mu=as.vector(exp(a %*% mvrnorm(1,coef(modfit1),vcov(modfit1)))),
         phi=modfit1$phi))*newdat$Effort
  }
  if(modtype=="TMBnbinom1") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$fit,
               TotalVar=.data$se.fit^2+.data$fit+.data$fit*sigma(modfit1))
      sim = replicate(nsim,simulateNegBin1Draw(modfit1,nObs,a,newdat$Effort))
  }
  if(modtype=="TMBnbinom2") {
       allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$fit,
               TotalVar=.data$se.fit^2+.data$fit+.data$fit^2/sigma(modfit1))
      sim = replicate(nsim,rnbinom(nObs,mu=exp(a %*% mvrnorm(1,fixef(modfit1)[[1]],
        vcov(modfit1)[[1]]))*newdat$Effort, size=sigma(modfit1)))
  }
  if(modtype=="TMBtweedie") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$Effort*.data$fit,
               TotalVar=.data$Effort^2*(.data$se.fit^2+sigma(modfit1)*.data$fit^(glmmTMB:::.tweedie_power(modfit1))))
       sim=replicate(nsim, simulateTMBTweedieDraw(modfit1,nObs,a,newdat$Effort) )
  }
    if(includeObsCatch & modtype!="Binomial") {
      d=match(allpred$matchColumn,obsdatval$matchColumn)
      allpred$Total[!is.na(d)]= allpred$Total[!is.na(d)] + obsdatval$Catch[d[!is.na(d)]]
      sim[!is.na(d),]= sim[!is.na(d),] + obsdatval$Catch[d[!is.na(d)]]
    }
    if(includeObsCatch & modtype=="Binomial") {
      d=match(allpred$matchColumn,obsdatval$matchColumn)
      allpred$Total[!is.na(d)]= obsdatval$pres[d[!is.na(d)]]
      sim[!is.na(d),]= obsdatval$presd[!is.na(d)]
    }
  stratatotal<-allpred %>%
      group_by_at(all_of(requiredVarNames)) %>%
      summarize(Total=sum(.data$Total,na.rm=TRUE))
  yeartotal<-allpred%>% group_by(.data$Year) %>%
      summarize(Total=sum(.data$Total,na.rm=TRUE))
  stratapred<-cbind(newdat,sim) %>%
       group_by_at(all_of(requiredVarNames)) %>%
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
  stratapred$Total=stratatotal$Total
  yearpred<-cbind(newdat,sim) %>%
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
  yearpred$Total<-yeartotal$Total
  if(is.na(max(yearpred$Total.cv)) | max(yearpred$Total.cv,na.rm=TRUE)>10) {
       print(paste(common[run],modtype," CV >10 or NA variance"))
       returnval=NULL
  }  else  {     returnval=yearpred  }
  if(printOutput) {
       write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modtype,"StratumSummary.csv"))
       write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modtype,"AnnualSummary.csv"))
  }
} else  {
       print(paste(common[run],modtype," CV >10 or NA variance"))
       returnval=NULL
}
  returnval
}

#' Generate predictions without estimating variance
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
#' @importFrom stats delete.response terms qnorm
#' @keywords internal
makePredictionsNoVar<-function(modfit1, modfit2=NULL, modtype, newdat, obsdatval=NULL, nsims, includeObsCatch, requiredVarNames, printOutput=TRUE, catchType, common, dirname, run) {
  if(includeObsCatch)    newdat$Effort=newdat$unsampledEffort/newdat$SampleUnits else
    newdat$Effort=newdat$Effort/newdat$SampleUnits
  newdat=uncount(newdat,.data$SampleUnits)
  getse=ifelse(modtype %in% c("Lognormal","Delta-Lognormal"),TRUE,FALSE)
  nObs=dim(newdat)[1]
  if(!is.null(modfit1)) {
    if(modtype=="Tweedie")     response1<-data.frame(cplm::predict(modfit1,newdata=newdat,type="response")) else
     response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=getse))
    if(dim(response1)[2]==1)     names(response1)="fit"
    if(!is.null(modfit2))  {
      response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=getse,type="response"))
      if(dim(response2)[2]==1) names(response2)[1]<-"fit"
      names(response2)=paste0(names(response2),"2")
    }
    if(modtype== "Binomial") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=.data$fit)
    }
    if(modtype== "Normal") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=.data$fit*.data$Effort)
    }
    if(modtype== "Lognormal") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=(lnorm.mean(.data$fit,sqrt(.data$se.fit^2+sigma(modfit1)^2))-0.1)*.data$Effort)
    }
    if(modtype== "Gamma") {
      allpred<-cbind(newdat,response1) %>%
        mutate(Total=(.data$fit-0.1)*.data$Effort)
    }
    if(modtype == "Delta-Lognormal" ){
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue=lnorm.mean(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2))) %>%
        mutate(Total=.data$Effort*.data$fit*.data$pos.cpue)
    }
    if(modtype == "Delta-Gamma"){
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(Total=.data$Effort*.data$fit*.data$fit2)
    }
    if(modtype =="NegBin") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$fit)
    }
    if(modtype =="Tweedie") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Total=.data$Effort*.data$fit)
    }
    if(modtype == "TMBnbinom1") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$fit)
    }
    if(modtype == "TMBnbinom2") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$fit)
    }
    if(modtype == "TMBtweedie") {
      allpred<-cbind(newdat,response1)  %>%
        mutate(Total=.data$Effort*.data$fit)
    }
    if(includeObsCatch & modtype!="Binomial") {
      a=match(allpred$matchColumn,obsdatval$matchColumn)
      allpred$Total[!is.na(a)]= allpred$Total[!is.na(a)] + obsdatval$Catch[a[!is.na(a)]]
    }
    if(includeObsCatch & modtype=="Binomial") {
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
        write.csv(stratapred,paste0(dirname[[run]],common[run],catchType[run],modtype,"StratumSummary.csv"))
        write.csv(yearpred,paste0(dirname[[run]],common[run],catchType[run],modtype,"AnnualSummary.csv"))
      }
  } else {
    returnval=NULL
  }
  returnval
}

#' Function to get an abundance index with SE, Does not yet have year interactions as random effects
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
#' @keywords internal
makeIndexVar<-function(modfit1, modfit2=NULL, modType, newdat, nsims, printOutput=FALSE, catchType = NULL, common = NULL, dirname = NULL, run = NULL) {
  returnval=NULL
  if(!is.null(modfit1)) {
    if(modType=="Tweedie")    response1<-data.frame(cplm::predict(modfit1,newdata=newdat,type="response",se.fit=TRUE)) else
    response1<-data.frame(predict(modfit1,newdata=newdat,type="response",se.fit=TRUE))
    if(dim(response1)[2]==1) {
      names(response1)="fit"
      if(modType=="Tweedie")
        response1$se.fit=getSimSE(modfit1, newdat, transFunc="exp", offsetval=NULL, nsim=nsims) else
          response1$se.fit=rep(NA,dim(response1)[2])
    }
    if(!is.null(modfit2))  {
      response2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
      names(response2)=paste0(names(response2),"2")
    }
    if(modType == "Delta-Lognormal" ){
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue=lnorm.mean(.data$fit2,.data$se.fit2),
               pos.cpue.se=lnorm.se(.data$fit2,.data$se.fit2),
               prob.se=.data$se.fit) %>%
        mutate(Index=.data$fit*.data$pos.cpue,
               SE=lo.se(.data$fit,.data$prob.se,.data$pos.cpue,.data$pos.cpue.se))
    }
    if(modType == "Delta-Gamma"){
      allpred<-cbind(newdat,response1,response2) %>%
        mutate(pos.cpue.se=.data$se.fit2,
               prob.se=.data$se.fit) %>%
        mutate(Index=.data$fit*.data$fit2,
               SE=lo.se(.data$fit,.data$prob.se,.data$fit2,.data$pos.cpue.se))
    }
    if(modType %in% c("Binomial","NegBin","Tweedie","TMBnbinom1","TMBnbinom2","Normal","TMBtweedie")) {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Index=.data$fit, SE=.data$se.fit)
    }
    if(modType =="Gamma") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Index=.data$fit-0.1, SE=.data$se.fit)
    }
    if(modType =="Lognormal") {
      allpred<-cbind(newdat,response1)   %>%
        mutate(Index=lnorm.mean(.data$fit,.data$se.fit)-0.1,
               SE=lnorm.se(.data$fit,.data$se.fit))
    }
    allpred=allpred %>% mutate(ymin=.data$Index-.data$SE,ymax=.data$Index+.data$SE)  %>%
      mutate(ymin=ifelse(.data$ymin<0,0,.data$ymin))
    returnval=allpred
    if(printOutput) {
      write.csv(allpred,paste0(dirname[[run]],common[run],catchType[run],modType,"Index.csv"))
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
                                            fittedPredictedResponse = predict(modfit1,type="response")))
        if(class(simulationOutput)[1]=="try-error") {
          simvals=simulateTweedie(modfit1,nsim*4)
          simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                              fittedPredictedResponse = predict(modfit1,type="response")))
        }
        if(class(simulationOutput)[1]=="try-error") {
          simvals=simulateTweedie(modfit1,nsim*10)
          simulationOutput = try(createDHARMa(simulatedResponse = simvals, observedResponse =modfit1$y ,
                                              fittedPredictedResponse = predict(modfit1,type="response")))
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
        g4<-ggplot(df1,aes(x=.data$Rank.Predictor,y=.data$Residual))+
          geom_point()+xlab("Model predictions (rank transformed)")+
          ylab("DHARMa scaled residuals")+ggtitle("d. Scaled residual vs. predicted")+
          geom_hline(aes(yintercept=0.5),lty=2)+geom_hline(aes(yintercept=0.75),lty=2)+geom_hline(aes(yintercept=0.25),lty=2)+
          geom_quantile(method = "rqss",col="red", formula=y ~ qss(x, lambda = 2))
        grid.arrange(g1,g2,g3,g4,ncol=2)
        test1=testUniformity(simulationOutput,plot=FALSE)
        test2=testDispersion(simulationOutput,plot=FALSE)
        test3=testZeroInflation(simulationOutput,plot=FALSE)
        test4=testOutliers(simulationOutput,plot=FALSE)
        returnval=c(test1$statistic,
                    test1$p.value,
                    test2$statistic,
                    test2$p.value,
                    test3$statistic,
                    test3$p.value,
                    test4$statistic,
                    test4$p.value)
        names(returnval)=c("KS.D","KS.p","Dispersion.ratio","Dispersion.p","ZeroInf.ratio","ZeroInf.p","Outlier","Outlier.p")
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
  if(modType=="Normal") {
    obsdatval$y=obsdatval$cpue
    modfit1<-try(lm(formula1,data=obsdatval))
  }
  if(modType=="Lognormal") {
    obsdatval$y=log(obsdatval$cpue+0.1)
    modfit1<-try(lm(formula1,data=obsdatval))
  }
  if(modType=="Delta-Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(lm(formula1,data=obsdatval))
  }
  if(modType=="Delta-Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit1=try(glm(formula1,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="NegBin") {
    obsdatval$y=round(obsdatval$Catch)
    modfit1=try(glm.nb(formula1,data=obsdatval,control=glm.control(epsilon=1E-6,maxit=45),na.action=na.fail))
  }
  if(modType=="Tweedie") {
    obsdatval$y=obsdatval$cpue
    modfit1=try(cpglm(formula1,data=obsdatval))
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
    if(modtype=="Tweedie")    predval1<-try(data.frame(cplm::predict(modfit1,newdata=newdat,type="response"))) else
      predval1<-try(data.frame(predict(modfit1,newdata=newdat,se.fit=TRUE,type="response")))
    if(class(predval1)[[1]]!="try-error") {
      if(!is.null(modfit2))  {
        predval2<-data.frame(predict(modfit2,newdata=newdat,se.fit=TRUE,type="response"))
        names(predval2)=paste0(names(predval2),"2")
      }
      if(modType=="Delta-Lognormal") {
        allpred<-cbind(newdat,predval1,predval2)  %>%
          mutate(est.cpue=.data$fit*lnorm.mean(.data$fit2,sqrt(.data$se.fit2^2+sigma(modfit2)^2)))
      }
      if(modType=="Delta-Gamma") {
        allpred<-cbind(newdat,predval1,predval2)   %>%
          mutate(est.cpue=.data$fit*.data$fit2)
      }
      if(modType %in% c("NegBin","TMBnbinom1","TMBnbinom2")) {
        allpred<-cbind(newdat,predval1)   %>%
          mutate(est.cpue=.data$fit/.data$Effort)
      }
      if(modType %in% c("TMBtweedie","Normal")){
        allpred<-cbind(newdat,predval1)   %>%
          mutate(est.cpue=.data$fit)
      }
      if(modType =="Gamma") {
        allpred<-cbind(newdat,predval1)   %>%
          mutate(est.cpue=.data$fit-0.1)
      }
      if(modType =="Tweedie") {
        allpred<-data.frame(est.cpue=predict(modfit1,newdata=newdat,type="response"))
      }
      if(modType =="Lognormal"){
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
  if(class(modfit)[1]=="cpglm") vcovval=modfit1$vcov else
      vcovval=vcov(modfit)
  if(length(coef(modfit))==dim(vcovval)[1]) {
    b=t(mvrnorm(nsim,coef(modfit),vcov(modfit)))
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

#' simulateNegBin1Draw
#'
#' @param modfit Value
#' @param nObs Value
#' @param b Value
#' @param Effort Value
#' @importFrom MASS mvrnorm
#' @importFrom glmmTMB fixef
#' @importFrom stats sigma rnbinom vcov
#' @keywords internal
simulateNegBin1Draw<-function(modfit,nObs,b,Effort) {
  muval<-exp(b %*% mvrnorm(1,fixef(modfit)[[1]],vcov(modfit)[[1]]))*Effort
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
#' @importFrom MASS mvrnorm gamma.shape
#' @importFrom glmmTMB fixef
#' @importFrom stats coef vcov rgamma
#' @keywords internal
simulateGammaDraw<-function(modfit,nObs,b) {
  muval<-exp(b %*% mvrnorm(1,coef(modfit),vcov(modfit)))
  shapeval<-gamma.shape(modfit)[[1]]
  scaleval<-muval/shapeval
  rgamma(nObs,shape=shapeval,scale=scaleval)
}

#' simulateTMBTweedieDraw
#'
#' @param modfit Value
#' @param nObs Value
#' @param b Value
#' @param Effort Value
#' @importFrom MASS mvrnorm
#' @importFrom stats vcov rgamma sigma
#' @importFrom tweedie rtweedie
#' @keywords internal
simulateTMBTweedieDraw<-function(modfit,nObs,b,Effort) {
  muval<-as.vector(exp(b %*% mvrnorm(1,fixef(modfit)[[1]],vcov(modfit)[[1]])))
  if(all(muval>0)) {
    simval<-rtweedie(nObs,power=glmmTMB:::.tweedie_power(modfit),
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
  cor12=cor(x1[!is.na(x1) &!is.na(x2)],x2[!is.na(x1) &!is.na(x2)])
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
plotSums<-function(yearpred,modType,fileName, subtext="", allVarNames, startYear, common, run, catchType, catchUnit) {
  if(is.numeric(yearpred$Year) & "Year" %in% allVarNames)
    yearpred$Year[yearpred$Source!="Ratio"]=yearpred$Year[yearpred$Source!="Ratio"]+startYear
  if(!is.null(yearpred)) {
    if(modType=="Binomial") ytitle=paste0(common[run]," ","predicted total positive trips") else
      ytitle=paste0("Total",common[run]," ",catchType[run]," (",catchUnit[run],")")
    yearpred<-yearpred %>%
      mutate(Year=as.numeric(as.character(.data$Year)),ymin=.data$Total-.data$Total.se,ymax=.data$Total+.data$Total.se) %>%
      mutate(ymin=ifelse(.data$ymin>0,.data$ymin,0))
    if(modType=="All") {
      g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI, fill=.data$Source))+
        geom_line(aes(col=.data$Source))+ geom_ribbon(alpha=0.3)+xlab("Year")+
        #        geom_line(aes(y=Total.mean,col=Source),lty=2,lwd=2)+
        ylab(ytitle)
    } else {
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
    if(modType=="Binomial") ytitle=paste0(common[run]," ","Positive trip index") else
      ytitle=paste0("Index ", common[run]," ",catchType[run]," (",catchUnit[run],")")
    if(modType %in% c("Delta-Lognormal","Delta-Gamma")) modType=paste("Delta",modType)
    yearpred<-yearpred %>% mutate(Year=as.numeric(as.character(.data$Year)),ymin=.data$Index-.data$SE,ymax=.data$Index+.data$SE) %>%
      mutate(ymin=ifelse(.data$ymin>0,.data$ymin,0))
    if(modType=="All") {
      g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Index,ymin=.data$ymin,ymax=.data$ymax,fill=.data$Source))+
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
plotSumsValidate<-function(yearpred,trueval,fileName,colName, allVarNames, startYear, common, run, catchType, catchUnit) {
  if(is.numeric(yearpred$Year)& "Year" %in% allVarNames) yearpred$Year[yearpred$Source!="Ratio"]=yearpred$Year[yearpred$Source!="Ratio"]+startYear
  yearpred<-yearpred %>%
    mutate(Year=as.numeric(as.character(.data$Year)),ymin=.data$Total-.data$Total.se,ymax=.data$Total+.data$Total.se) %>%
    mutate(ymin=ifelse(.data$ymin>0,.data$ymin,0))
  trueval<-trueval %>% rename(Total=!!colName) %>%
    mutate(ymin=NA,ymax=NA,Source="Validation",
           Total.mean=NA,TotalLCI=NA,TotalUCI=NA)
  yearpred<-bind_rows(yearpred[,c("Year","Total","Total.mean","TotalLCI","TotalUCI","ymin","ymax","Source")],
                      trueval[,c("Year","Total","Total.mean","TotalLCI","TotalUCI","ymin","ymax","Source")])
  if(all(is.na(yearpred$Total.mean)))
    g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI,fill=.data$Source))+
    geom_line(aes(color=.data$Source))+ geom_ribbon(alpha=0.3)+
    xlab("Year")+
    ylab(paste0(common[run]," ",catchType[run]," (",catchUnit[run],")"))+
    geom_point(data=yearpred[yearpred$Source=="Validation",],aes(x=.data$Year,y=.data$Total,color=.data$Source),size=2) else
      g<-ggplot(yearpred,aes(x=.data$Year,y=.data$Total,ymin=.data$TotalLCI,ymax=.data$TotalUCI,fill=.data$Source))+
    geom_line(aes(color=.data$Source))+ geom_ribbon(alpha=0.3)+
    geom_line(aes(y=.data$Total.mean,color=.data$Source),lty=2)+
    xlab("Year")+
    ylab(paste0(common[run]," ",catchType[run]," (",catchUnit[run],")"))+
    geom_point(data=yearpred[yearpred$Source=="Validation",],aes(x=.data$Year,y=.data$Total,color=.data$Source),size=2)
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
  if(modType=="Delta-Lognormal") {
    obsdatval$y=obsdatval$log.cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit2=try(lm(formula2,data=obsdatval))
  }
  if(modType=="Delta-Gamma") {
    obsdatval$y=obsdatval$cpue
    obsdatval=obsdatval[obsdatval$cpue>0,]
    modfit2=try(glm(formula2,data=obsdatval,family=Gamma(link="log")))
  }
  if(modType=="NegBin") {
    obsdatval$y=round(obsdatval$Catch)
    modfit1=try(glm.nb(formula3,data=obsdatval,control=glm.control(epsilon=1E-6,maxit=45),na.action=na.fail))
  }
  if(modType=="Tweedie") {
    obsdatval$y=obsdatval$cpue
    modfit1=try(cpglm(formula2,data=obsdatval))
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
  if(class(modfit1)[1]=="try-error") modfit1=NULL
  if(class(modfit2)[1]=="try-error") modfit2=NULL
  if(!is.null(modfit1)) {
    if(modType %in% c("Binomial","Delta-Lognormal","Delta-Gamma"))  #for delta models write binomial anova
      write.csv(anova(modfit1,test="Chi"),file=paste0(outputDir,"/BinomialAnova.csv"))
    if(modType %in% c("NegBin")) anova1=anova(modfit1,test="Chi")
    if(modType %in% c("Tweedie","TMBtweedie","TMBnbinom1","TMBnbinom2")) anova1=NULL
    if(modType %in% c("Delta-Lognormal","Delta-Gamma")) anova1=anova(modfit2,test="F")
    if(!is.null(anova1)) {
      write.csv(anova1,paste0(outputDir,"/anova",modType,".csv"))
    }
  }
  list(modfit1=modfit1,modfit2=modfit2)
}


#' Function to make data summarizes including ratio estimate at
#'
#'stratification defined by strataVars
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


