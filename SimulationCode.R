###############################################
#Setting up simulation for continuous biomarker
###############################################
#Adapting simulation from Riz 2017 Biom J paper
#Coded LM models for prediction 
#Next: Code PE for LM model 
#Next: Code Copula models for prediction
#Next: Code PE for Copula models
#Next: Code PE for Joint models

rm(list=ls())

#Load in required packages
library(JM)
library(survival)
library(MASS)
library(dynpred)
library(reshape2)
library(timeROC)
library(pec)
library(pbivnorm)
library(lmvar)

#Set working directory to where "PredictionFunctions.R" is located 
# setwd()

seed<-17 #Change for each simulation 
set.seed(seed)

#parameters to change: 
sigma.y <- 0.6 # measurement error standard deviation #0.6, 1.2
alpha <- 1.5 # association parameter #0.5, 1.5

int.base<-ifelse(alpha==0.5,-4,-6) #-4, -6
insp.rate<-0.5 #1 or 0.5 for more sparse observations
thor<-3

N <- 1000 # number of subjects
K <- 15  # number of planned repeated measurements per subject, per outcome
t.max <- 15 # maximum follow-up time
cens_horiz <-15 

################################################

# parameters for the linear mixed effects model
betas <- c("(Intercept)"=-3,
           "time"=1,
           "Group" = -0.8, 
           "Group:time" = 0.5)


# parameters for the survival model
gammas <- c("(Intercept)" = int.base, "Group" = 0.5)  # coefficients for baseline covariates

phi <- 1.4 # shape for the Weibull baseline hazard

D <- matrix(c(1, 0.5, 
              0.5, 1), 2, 2)
D <- (D + t(D)) / 2

################################################

# at which time points longitudinal measurements are supposed to be taken (Fixed every year)
# times <- c(replicate(N, c(0, sort(runif(K-1, 0, t.max)))))
# times <- c(replicate(N, 0:(K-1)))
#Inspections are made using an exponential distribution with rate=insp.rate
times <- c(replicate(N,cumsum(c(0,rexp(n=K-1,rate=insp.rate)))))
times_dat<-data.frame(id=rep(1:N,each=K),time=times)

group <- rbinom(N, 1, 0.5) # group indicator (X)
group_dat<-data.frame(id=1:N,group)

DF<-merge(times_dat,group_dat,by="id",all.x=TRUE)

# design matrices for the longitudinal measurement model
X <- model.matrix(~ time + group + group:time, data = DF)
Z <- model.matrix(~ time, data = DF)

# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Group" = group)

################################################
max.FUtime <- max(times) + 2 * IQR(times)

#simulate random effects
b <- mvrnorm(N, rep(0, nrow(D)), D)

# simulate longitudinal responses
id <- rep(1:N, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ]))
y <- rnorm(N * K, eta.y, sigma.y)


# simulate event times
eta.t <- as.vector(W %*% gammas)
invS <- function (t, u, i) {
  h <- function (s) {
    group <- group[i]
    XX <- cbind(1, s, group, group*s)
    ZZ <- cbind(1, s)
    f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha)
  }
  integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(N)
trueTimes <- numeric(N)
for (i in 1:N) {
  Root <- try(uniroot(invS, interval = c(1e-05, max.FUtime), u = u[i],
                      i = i,extendInt="upX")$root, TRUE)
  #"Cure" patients have event time set to Infinity 
  trueTimes[i] <- ifelse(inherits(Root, "try-error"),Inf,Root)
}

na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
group <- group[na.ind]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

#Simulate censoring times from uniform distribution 
Ctimes<-runif(n,0,cens_horiz)
Time <- pmin(trueTimes, Ctimes,rep(cens_horiz,n))
event <- ifelse(trueTimes<=Time,1,0)
prop.table(table(event))*100

Time_dat<-data.frame(id=rep(1:N,each=K),Time=rep(Time,each=K))

################################################
# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= Time_dat$Time
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
dat.id<-dat[!duplicated(dat$id),]

#Add in the baseline marker value 
dat<-merge(x=dat,y=dat.id[,c("id","y")],by="id",suffixes=c("","0"))

#Select 500 for training set (500 for test)
set <- sort(sample(unique(id), 500))
train_data <-subset(dat,!dat$id%in%set) #dat[!dat$id %in% set, ]
test_data <-subset(dat,dat$id%in%set) # dat[dat$id %in% set, ]

trueValues <- list(betas = betas, phi = phi, gammas = gammas, alpha = alpha,
                   b_test = b) 

train_data.id<-train_data[!duplicated(train_data$id),]
test_data.id<-test_data[!duplicated(test_data$id),]

###############################################################################
###############################################################################
### Joint Models
###############################################################################
###############################################################################
# True joint model
lmeFit<-lme(y~group+time+group:time,data=train_data,
            random=~time|id)
survFit <- coxph(Surv(Time, event) ~ group, data = train_data.id,
                 x = TRUE)
jointFit1 <- jointModel(lmeFit, survFit, timeVar = "time")

# Misspecified joint model 
lmeFit2<-lme(y~group+I(time^3)+group:I(time^3),data=train_data,
             random=~time|id)
survFit2 <- coxph(Surv(Time, event) ~ group, data = train_data.id,
                  x = TRUE)
jointFit2 <- try(jointModel(lmeFit2, survFit2, timeVar = "time"),TRUE)

###############################################################################
###############################################################################
### Landmark Models
###############################################################################
###############################################################################
LMdata <- NULL
LMs <- seq(0,7,by=1)

LMdata<-cutLM(data=train_data.id,outcome=list(time="Time",status="event"),
              LM=0,horizon=thor,covs=list(fixed=c("group"),varying="y"),
              format="long",id="id",rtime="time")

LMdata$y<-subset(train_data,train_data$time==0)$y
LMdata$time<-0

for (i in 2:length(LMs)) {
  LMdata <- rbind(LMdata,cutLM(data=train_data,outcome=list(time="Time",status="event"),
                               LM=LMs[i],horizon=LMs[i]+thor,covs=list(fixed=c("group"),varying="y"),
                               format="long",id="id",rtime="time",right=FALSE))
}

LMdata$survstatus<-LMdata$event
LMdata$wsurvtime<-LMdata$Time

tt<-sort(unique(LMdata$wsurvtime[LMdata$survstatus==1]))
LMdata$Tstart<-LMdata$LM
LMdata2<-survSplit(Surv(Tstart,wsurvtime,survstatus)~.,data=LMdata,cut=tt,end="wsurvtime",start="Tstart",event="survstatus")

LMdata2$y_nph <- LMdata2$y*(LMdata2$wsurvtime-LMdata2$LM)
LMdata2$y_nph2 <- LMdata2$y*(LMdata2$wsurvtime-LMdata2$LM)^2

LMdata2$y_tau<-LMdata2$y*LMdata2$LM
LMdata2$y_tau2<-LMdata2$y*LMdata2$LM^2

# ipl*
g1 <- function(t) (t)
g2 <- function(t) (t)^2

LMdata2$LM1 <- g1(LMdata2$LM)
LMdata2$LM2 <- g2(LMdata2$LM)

LMdata2$group_int<-LMdata2$y*LMdata2$group
LMdata$group_int<-LMdata$y*LMdata$group

####################################
#LM models (stratified) from superdataset
LMsupercox0 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group +  y + strata(LM) + cluster(id), data=LMdata, method="breslow")
bh_supercox0<-basehaz(LMsupercox0,centered=FALSE)
LMsupercox1 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 + strata(LM) + cluster(id),
                     data=LMdata2, method="breslow")
bh_supercox1<-basehaz(LMsupercox1,centered=FALSE)


LMsupercox0_int <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y +
                           group_int +
                           strata(LM) + cluster(id),
                         data=LMdata, method="breslow")
bh_supercox0_int<-basehaz(LMsupercox0_int,centered=FALSE)

LMsupercox1_int <-coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 +
                          group_int +
                          strata(LM) + cluster(id),
                        data=LMdata2, method="breslow")
bh_supercox1_int<-basehaz(LMsupercox1_int,centered=FALSE)

####################################
#Restructure LMdataset at longitudinal 
#Create LMdata set with administrative censoring at LM+w
long_data<-train_data
long_data$wsurvtime<-pmin(long_data$time+thor,long_data$Time)
long_data$survstatus<-ifelse(long_data$Time<=long_data$wsurvtime,1,0)
LMdata<-long_data
tt<-sort(unique(LMdata$wsurvtime[LMdata$survstatus==1]))
LMdata$LM<-LMdata$time
LMdata$Tstart<-LMdata$time
LMdata2<-survSplit(Surv(Tstart,wsurvtime,survstatus)~.,data=LMdata,cut=tt,end="wsurvtime",start="Tstart",event="survstatus")

LMdata2$y_nph <- LMdata2$y*(LMdata2$wsurvtime-LMdata2$LM)
LMdata2$y_nph2 <- LMdata2$y*(LMdata2$wsurvtime-LMdata2$LM)^2

LMdata2$y_tau<-LMdata2$y*LMdata2$LM
LMdata2$y_tau2<-LMdata2$y*LMdata2$LM^2


# ipl*
g1 <- function(t) (t)
g2 <- function(t) (t)^2

LMdata2$LM1 <- g1(LMdata2$LM)
LMdata2$LM2 <- g2(LMdata2$LM)

LMdata2$group_int<-LMdata2$y*LMdata2$group
LMdata$group_int<-LMdata$y*LMdata$group
####################################
#LM models with interactions 
LMsupercox2_int <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau +y_tau2 +
                           group_int + 
                           LM1 + LM2 + cluster(id), data=LMdata2, method="breslow")
bh_supercox2_int<-basehaz(LMsupercox2_int,centered=FALSE)

LMsupercox3_int <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_nph + y_nph2 + 
                           group_int + 
                           LM1 + LM2 + cluster(id), data=LMdata2, method="breslow")
bh_supercox3_int<-basehaz(LMsupercox3_int,centered=FALSE)

LMsupercox4_int <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 + y_nph + y_nph2 + 
                           group_int + 
                           LM1 + LM2 + cluster(id), 
                         data=LMdata2, method="breslow")
bh_supercox4_int<-basehaz(LMsupercox4_int,centered=FALSE)

LMsupercox2 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau +y_tau2 +
                       LM1 + LM2 + cluster(id), data=LMdata2, method="breslow")
bh_supercox2<-basehaz(LMsupercox2,centered=FALSE)
LMsupercox3 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_nph + y_nph2 + LM1 + LM2 + cluster(id),
                     data=LMdata2, method="breslow")
bh_supercox3<-basehaz(LMsupercox3,centered=FALSE)
LMsupercox4 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 + y_nph + y_nph2 +
                       LM1 + LM2 + cluster(id),
                     data=LMdata2, method="breslow")
bh_supercox4<-basehaz(LMsupercox4,centered=FALSE)

##########################

###############################################################################
###############################################################################
### Landmark Models
###############################################################################
###############################################################################
#Step 1: Model for T 
################################################
#Model for T 
cox_T<-coxph(Surv(Time,event)~group,data=train_data.id)
cox_T_base<-coxph(Surv(Time,event)~group+y0,data=train_data.id)
#Test nonPH
zp<-cox.zph(cox_T)

#Non-proportional hazards models
dtimes<-sort(unique(with(train_data.id,Time[event==1])))
tdata<-survSplit(Surv(Time,event)~.,data=train_data.id,cut=dtimes)
tdata$y0_lt<-log(tdata$Time+1)*tdata$y0
tdata$y0_t<-tdata$Time*tdata$y0
tdata$y0_t2<-(tdata$Time^2)*tdata$y0

cox_T_nph<-coxph(Surv(tstart,Time,event)~group+y0+y0_t+y0*group,data=tdata)
cox_T_nph_sq<-coxph(Surv(tstart,Time,event)~group+y0+y0_t+y0_t2+y0*group,data=tdata)
cox_T_nph_logT<-coxph(Surv(tstart,Time,event)~group+y0+y0_lt+y0*group,data=tdata)

#Parametric Weibull hazard 
weib_T<-survreg(Surv(Time,event)~group,data=train_data.id,dist="weibull")
#   survreg's scale  =    1/(rweibull shape)
#   survreg's intercept = log(rweibull scale)
shape.weib<-1/weib_T$scale
scale.weib<-exp(weib_T$coef[1])
coef.weib<- -weib_T$coef[-1]/weib_T$scale

#Function that returns F values for all cox models considered
#NOTE: Takes long time to run with the inclusion of non proprotional hazards models (choose models and modify code accordingly)
F_T_func_cond<-function(dat)
{
  xdata1<-as.numeric(dat[1])
  tt<-as.numeric(dat[2])
  tt_LM<-as.numeric(dat[3])

  bh<-list(basehaz(cox_T,FALSE)
           # ,
           # basehaz(cox_T_nph,FALSE),
           # basehaz(cox_T_nph_sq,FALSE),
           # basehaz(cox_T_nph_logT,FALSE)
  )
  
  bet<-list(cox_T$coef
            # ,
            # cox_T_nph$coef,
            # cox_T_nph_sq$coef,
            # cox_T_nph_logT$coef
            )
  names(bet)<-c("cox_T"
                # ,
                # "cox_T_nph","cox_T_nph_sq","cox_T_nph_logT"
                )
  
  sfi<-bh
  sfi<-lapply(sfi,function(x) data.frame(x,haz0=diff(c(0,x$hazard))))
  
  sfi[[1]]["haz"]<-sfi[[1]]$haz0*as.numeric(exp(bet[["cox_T"]][1]*xdata1))
  # sfi[[2]]["haz"]<-sfi[[2]]$haz0*as.numeric(exp(bet[["cox_T_nph"]][1]*xdata1+
  #                                                 bet[["cox_T_nph"]][2]*y0+
  #                                                 bet[["cox_T_nph"]][3]*y0*sfi[[2]]$time+
  #                                                 bet[["cox_T_nph"]][4]*y0*xdata1))
  # sfi[[3]]["haz"]<-sfi[[3]]$haz0*as.numeric(exp(bet[["cox_T_nph_sq"]][1]*xdata1+
  #                                                 bet[["cox_T_nph_sq"]][2]*y0+
  #                                                 bet[["cox_T_nph_sq"]][3]*y0*sfi[[3]]$time+
  #                                                 bet[["cox_T_nph_sq"]][4]*y0*(sfi[[3]]$time)^2+
  #                                                 bet[["cox_T_nph_sq"]][5]*y0*xdata1))
  # sfi[[4]]["haz"]<-sfi[[4]]$haz0*as.numeric(exp(bet[["cox_T_nph_logT"]][1]*xdata1+
  #                                                 bet[["cox_T_nph_logT"]][2]*y0+
  #                                                 bet[["cox_T_nph_logT"]][3]*y0*log(sfi[[4]]$time+1)+
  #                                                 bet[["cox_T_nph_logT"]][4]*y0*xdata1))
  
  sfi<-lapply(sfi,function(x) data.frame(x,Haz=cumsum(x$haz)))
  tmp<-lapply(sfi,function(x) evalstep(x$time,x$Haz,c(tt,tt_LM),subst=0))
  Fw<-lapply(tmp,function(x) ((1-exp(-x[1]))-(1-exp(-x[2])))/exp(-x[2]))  #-exp(-x[1])/exp(-x[2])+1) #(-exp(-x[1])+exp(-x[2]))/(exp(-x[2])))
  
  Fw_weib_tt<-exp(-(tt/scale.weib)^shape.weib*exp(coef.weib*xdata1))
  Fw_weib_tt_LM<-exp(-(tt_LM/scale.weib)^shape.weib*exp(coef.weib*xdata1))
  Fw_weib<-((1-Fw_weib_tt)-(1-Fw_weib_tt_LM))/Fw_weib_tt_LM
  
  return(c(unlist(Fw),Fw_weib))
}

FT_cond_vec<-apply(LMdata[,c("group","Time","LM")],1,F_T_func_cond)
FT_cond_vec<-t(FT_cond_vec)
FT_names<-c("cox",
            # "nph","nphsq","nphlog",
            "weib") 
FT_cond_vec<-data.frame(FT_cond_vec)
colnames(FT_cond_vec)<-FT_names
#FT_cond_vec contains the conditional survival predictions for patients in the training data set for all of the survival models considered (different columns)

################################################
#Step 2: Model for Z (no longer need model for Z*!) 
################################################
F_Zstar<-function(zstar,mean_Zstar,sigma_Zstar,log.prob=FALSE)
{
  pnorm(zstar,mean=mean_Zstar,sd=sigma_Zstar,log.p=log.prob)
}

f_Zstar<-function(zstar,mean_Zstar,sigma_Zstar)
{
  dnorm(zstar,mean=mean_Zstar,sd=sigma_Zstar)
}

#Models for Z (with intercept)
const<-~1
simp<-~LM+group
int<-~LM*group
sq<-~LM+I(LM^2)+group
sqint<-~LM*group+I(LM^2)*group
sqpartint<-~LM*group+I(LM^2)
bsp<-~bs(LM, Boundary.knots=c(0,cens_horiz))+group
bspint<-~bs(LM, Boundary.knots=c(0,cens_horiz))*group

mods<-list(const,simp,int,sq,sqint,sqpartint,bsp,bspint)
mods_names<-c("const","simp","int","sq","sqint","sqpartint","bsp","bspint")
names(mods)<-mods_names

#Models for Z (without intercept)
const.red<-~1-1
simp.red<-~LM+group-1
int.red<-~LM*group-1
sq.red<-~LM+I(LM^2)+group-1
sqint.red<-~LM*group+I(LM^2)*group-1
sqpartint.red<-~LM*group+I(LM^2)-1
bsp.red<-~bs(LM, Boundary.knots=c(0,cens_horiz))+group-1
bspint.red<-~bs(LM, Boundary.knots=c(0,cens_horiz))*group-1

mods.red<-list(const.red,simp.red,int.red,sq.red,sqint.red,sqpartint.red,bsp.red,bspint.red)
names(mods.red)<-c("const","simp","int","sq","sqint","sqpartint","bsp","bspint")

#Chose models to fit for mu, sigma, and rho (the association function)
mean_names<-c("simp","int","sqint","bsp","bspint") #models fit for mean of Z 
sd_names<-c("const","simp","int","bsp","bspint") #models fit for sigma of Z
rho_names<-c("const","simp","int","bsp","bspint") #models fit for association function 

#Fitting model for Z allowing mead and variance to be the functions considered above
#NOTE: choose which models you want to fit to minimize running time 
for(mean_string in mean_names)
{
  for(sd_string in sd_names)
  {
    mod_mean_Zstar<-mods.red[[mean_string]]
    mod_sd_Zstar<-mods.red[[sd_string]]
    
    xmat.mean<-model.matrix(mod_mean_Zstar,data=LMdata)
    xmat.var<-model.matrix(mod_sd_Zstar,data=LMdata)
    
    fit<-lmvar(LMdata$y,xmat.mean,xmat.var)
    assign(paste0("mod_Zstar_full_mean",mean_string,"_sd",sd_string),fit)
    assign(paste0("mod_Zstar_mean",mean_string,"_sd",sd_string),coef(fit))
  }
}

##################################
#Joint distribution
##################################
FZstarT<-function(FT_cond,zstar,mean_Zstar,sigma_Zstar,rho) #,df)
{
  a<-qnorm(FT_cond,mean=0,sd=1)
  b<-qnorm(F_Zstar(zstar,mean_Zstar,sigma_Zstar),mean=0,sd=1)
  ret<-ifelse(a==-Inf,0,
              ifelse(a==Inf,F_Zstar(0,mean_Zstar,sigma_Zstar),
                     ifelse(b==Inf,FT_cond,pbivnorm(x=a,y=b,rho=rho))))
  
  #Uncomment for fitting Student's t copula 
  # a<-qt(FT_cond,df=df)
  # b<-qt(F_Zstar(zstar,mean_Zstar,sigma_Zstar),df=df)
  # ret<-ifelse(a==-Inf,0,
  #             ifelse(a==Inf,F_Zstar(0,mean_Zstar,sigma_Zstar),
  #                    ifelse(b==Inf,FT_cond,pmvt(lower=c(-Inf,-Inf),upper=c(a,b),corr=rho*matrix(c(1,0,0,1),nrow=2,df=df)))))
  return(ret)
}

##################################
#Maximize joint likelihood
##################################
#Association models
loglik_ZT<-function(g,FTcond,meanZ,sigmaZ,dat,y) #,v) #Uncomment for Student's t copula 
{
  eta<-coef_dat%*%g
  rho<-1-2/(1+exp(2*eta))
  
  q1<-qnorm(FTcond,mean=0,sd=1)
  q2<-qnorm(F_Zstar(y,meanZ,sigmaZ,log.prob=TRUE),mean=0,sd=1,log.p=TRUE)
  
  #Student's t copula
  # q1<-qt(FTcond,df=v)
  # q2<-qt(F_Zstar(dat$y,meanZ,sigmaZ,log.prob=TRUE),df=v,log.p=TRUE)
  
  #L1: P(T=t,Z=z)=f(t,z)
  L1<--1/2*log(1-rho^2)-(rho^2*(q1^2+q2^2)-2*rho*q1*q2)/(2*(1-rho^2))
  
  #Student's t copula 
  #L1<--(v+1)/2*log(1+1/(v*(1-rho^2))*(q1^2+q2^2-2*rho*q1*q2))-dt(q1,df=v,log=TRUE)-dt(q2,df=v,log=TRUE)-1/2*log(1-rho^2)
  
  L1_cont<-L1[which(dat$event==1)]
  
  #L2: P(T>t,Z=z)
  L2<-log(pnorm(-(q1-rho*q2)/sqrt(1-rho^2)))
  
  #Student's t copula
  #L2<-log(pt(-(q1-rho*q2)/sqrt((v+q2^2)*(1-rho^2)/(v+1)),df=v+1))
  
  L2_cont<-L2[which(dat$event==0)]
  
  ret<--(sum(L1_cont)+sum(L2_cont))
  
  return(ret)
}

#Perform maximization for the different copula models
#Error for one-dimensional optimization if "const" model considered 
for(rho_string in rho_names)
{
  for(mean_string in mean_names)
  {
    for(sd_string in sd_names)
    {
      for(FT_string in FT_names)
      {
        print(c(rho_string,mean_string,sd_string,FT_string))
        xmat.mean.pred<-model.matrix(mods.red[[mean_string]],LMdata)
        xmat.sd.pred<-model.matrix(mods.red[[sd_string]],LMdata)
        modZ_vals<-predict(get(paste0("mod_Zstar_full_mean",mean_string,"_sd",sd_string)),xmat.mean.pred,xmat.sd.pred)
        meanZ<-modZ_vals[,1]
        sdZ<-modZ_vals[,2]
        
        coef_dat<-model.matrix(mods[[rho_string]],data=LMdata)
        
        start_val<-rep(0,ncol(coef_dat))
        temp<-try(optim(start_val,loglik_ZT,FTcond=FT_cond_vec[,FT_string],
                        meanZ=meanZ,
                        sigmaZ=sdZ,
                        y=LMdata$y,
                        # v=df0, #Uncomment for Student's t copula 
                        dat=LMdata,
                        method="BFGS",control=list(maxit=5000,trace=FALSE)))
        
        if(!inherits(temp,"try-error"))
        {
          conv_val<-temp$value
          conv_tol<-1
          while(conv_tol>10e-5)
          {
            temp<-optim(temp$par,loglik_ZT,FTcond=FT_cond_vec[,FT_string],
                        meanZ=meanZ,
                        sigmaZ=sdZ,
                        y=LMdata$y,
                        # v=df0, #Uncomment for Student's t copula 
                        dat=LMdata,
                        method="Nelder-Mead",control=list(maxit=5000,trace=FALSE)) #,hessian=TRUE)
            conv_tol<-abs(temp$value-conv_val)
            conv_val<-temp$value
          }
          
          assign(paste0("MLE_rho",rho_string,"_mean",mean_string,"_sd",sd_string,"_",FT_string),temp$par) 
        } else {
          print(c(rho_string,mean_string,sd_string,FT_string,"ERROR"))
          assign(paste0("MLE_rho",rho_string,"_mean",mean_string,"_sd",sd_string,"_",FT_string),start_val) 
        }
      }
    }
  }
}


##############################################
#Evaluating predictions  
source("PredictionFunctions.R")

#Prediction horizon
w_predict<-thor

#Names of models to fit
models<-c("Null","Truth","LM0","LM1","LMInt0","LMInt1","JM","JM1_sim","JM2","JM2_sim","CC1","CW1")
nmodels<-length(models)

#Set vectors to be the copula models to be fit (mean, sd, rho)
mean_names_fit<-c("bspint")
sd_names_fit<-c("simp")
rho_names_fit<-c("int")

#Lanmdark times at which to evaluate predictions 
LMx<-seq(0,5,by=1)

BS_full<-NULL 

for(t0 in LMx)
{
  print(t0)
  
  pred1<-summary(survfit(Surv(Time,event)~1,data=train_data.id),time=t0)$surv
  pred2<-summary(survfit(Surv(Time,event)~1,data=train_data.id),time=t0+w_predict)$surv
  pred.full<-(pred1-pred2)/pred1
  
  #Subset of data with those who are still alive at time t0
  sub_dat<-subset(test_data,test_data$Time>t0)
  sub_dat.id<-subset(test_data.id,test_data.id$Time>t0)
  
  #select those with event time>LM and their data up to LM 
  sub<-subset(sub_dat,sub_dat$time<=t0) 
  sub<-sub[order(sub$id,sub$time),]
  #select last entry for LM models 
  sub_LM<-do.call("rbind", as.list(by(sub,sub$id,tail,n=1)))
  sub_LM$LM<-t0
  
  data.temp<-rep(pred.full,nrow(sub_dat.id))
  
  zdata=sub_LM$y
  xdata=sub_LM$group 
  # ydata=sub_LM$y0
  obs.time<-sub_LM$time
  data.temp<-cbind(data.temp,BSpredict_Truth(trueValues,sub_data=sub_dat.id,t0=t0,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM0(bh=bh_supercox0,bet=LMsupercox0$coef,zdata=zdata,t0=t0,xdata=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LM1(bh=bh_supercox1,bet=LMsupercox1$coef,zdata=zdata,t0=t0,xdata=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LMInt0(bh=bh_supercox0_int,bet=LMsupercox0_int$coef,zdata=zdata,t0=t0,xdata=xdata,w_predict))
  data.temp<-cbind(data.temp,BSpredict_LMInt1(bh=bh_supercox1_int,bet=LMsupercox1_int$coef,zdata=zdata,t0=t0,xdata=xdata,w_predict))
  if(inherits(jointFit1,"try-error"))
  {
    data.temp<-cbind(data.temp,NA,NA)
  } else {
    data.temp<-cbind(data.temp,BSpredict_JM(JM_mod=jointFit1,sub_data=sub,t0=t0,w_predict,sim=TRUE))
    data.temp<-cbind(data.temp,BSpredict_JM(JM_mod=jointFit1,sub_data=sub,t0=t0,w_predict,sim=FALSE)) #Takes long to run
  }
  if(inherits(jointFit2,"try-error"))
  {
    data.temp<-cbind(data.temp,NA,NA)
  } else {
    data.temp<-cbind(data.temp,BSpredict_JM(JM_mod=jointFit2,sub_data=sub,t0=t0,w_predict,sim=TRUE))
    data.temp<-cbind(data.temp,BSpredict_JM(JM_mod=jointFit2,sub_data=sub,t0=t0,w_predict,sim=FALSE)) #Takes long to run
  }
  F_T_dat<-apply(cbind(xdata,t0+w_predict,t0),1,F_T_func_cond)
  F_T_dat<-data.frame(t(F_T_dat))
  names(F_T_dat)<-FT_names
  
  for(mean_string in mean_names_fit) 
  {
    for(sd_string in sd_names_fit)
    {
      for(rho_string in rho_names_fit) 
      {
        data.temp<-cbind(data.temp,BSpredict_Copula(zdata=zdata,t0=t0,
                                                    # v=df0, #Uncomment if using Student's t copula 
                                                    xdata=xdata,
                                                    w_predict=w_predict,mean_string=mean_string,sd_string=sd_string,
                                                    F_T_dat=F_T_dat,FT_string=FT_names,rho_string=rho_string,obs.time=obs.time))
      }
    }
  }
  
  df.temp<-data.frame(sub_LM,data.temp)
  
  BS_full<-rbind(BS_full,df.temp)
}

names(BS_full)<-c(names(sub_LM),models)
#BS_full: contains the predictions for each of the models considered for each of the individuals still alive at that time in the 
#test dataset. 

#########################
#RMSE
#########################
RMSE_all<-sqrt(colMeans((BS_full-BS_full$Truth)^2))[models]

RMSE_X<-NULL
RMSE_LM<-NULL
RMSE_LM_X<-NULL
for(i in LMx)
{
  RMSE_sub<-subset(BS_full,BS_full$LM==i)
  RMSE_LM<-rbind(RMSE_LM,c(LM=i,sqrt(colMeans((RMSE_sub-RMSE_sub$Truth)^2))[models]))
  for(j in c(0,1))
  {
    RMSE_sub<-subset(BS_full,BS_full$group==j&BS_full$LM==i)
    RMSE_LM_X<-rbind(RMSE_LM_X,c(LM=i,group=j,sqrt(colMeans((RMSE_sub-RMSE_sub$Truth)^2))[models]))
    
    RMSE_sub<-subset(BS_full,BS_full$group==j)
    RMSE_X<-rbind(RMSE_X,c(group=j,sqrt(colMeans((RMSE_sub-RMSE_sub$Truth)^2))[models]))
  }
}
#RMSE_LM_X: RMSE by baseline covariate and landmark time 
#RMSE_X: RMSE by baseline covariate X 

#########################
#AUC and BS 
#########################
df_BS<-data.frame(matrix(nrow=length(LMx),ncol=nmodels))
df_AUC<-data.frame(matrix(nrow=length(LMx),ncol=nmodels))
for(j in 1:nmodels)
{
  if(all(is.na(BS_full[models[j]]))==TRUE)
  {
    df_BS[,j]<-NA
    df_AUC[,j]<-NA
  } else {
    pred.error<-PE(models[j],LMx,thor,BS_full)
    df_BS[,j]<-pred.error[[1]]
    df_AUC[,j]<-pred.error[[2]]
  }
}

names(df_BS)<-models
names(df_AUC)<-models
#df_AUC: AUC for each the models by landmark time 
#df_BS: Brier score for each of the models by landmark time 