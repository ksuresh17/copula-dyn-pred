###############################################
#Setting up simulation for continuous biomarker
###############################################
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

###############################################################################
###############################################################################
### Landmark Models
###############################################################################
###############################################################################
LMdata <- NULL
LMs <- seq(0,5,by=1)

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

#Extended
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
LMsupercox1 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 + strata(LM) + cluster(id),
                     data=LMdata2, method="breslow")
bh_supercox1<-basehaz(LMsupercox1,centered=FALSE)

LMsupercox1_int <-coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 +
                          group_int +
                          strata(LM) + cluster(id),
                        data=LMdata2, method="breslow")
bh_supercox1_int<-basehaz(LMsupercox1_int,centered=FALSE)

####################################
#Longitudinal landmark data set for supermodels
LMdata_long<-train_data
LMdata_long$wsurvtime<-pmin(LMdata_long$time+thor,LMdata_long$Time)
LMdata_long$survstatus<-ifelse(LMdata_long$Time<=LMdata_long$wsurvtime,1,0)

#Extended 
tt<-sort(unique(LMdata_long$wsurvtime[LMdata_long$survstatus==1]))
LMdata_long$LM<-LMdata_long$time
LMdata_long$Tstart<-LMdata_long$LM
LMdata2_long<-survSplit(Surv(Tstart,wsurvtime,survstatus)~.,data=LMdata_long,cut=tt,end="wsurvtime",start="Tstart",event="survstatus")

LMdata2_long$y_tau<-LMdata2_long$y*LMdata2_long$LM
LMdata2_long$y_tau2<-LMdata2_long$y*LMdata2_long$LM^2

#ipl*
LMdata2_long$LM1 <- g1(LMdata2_long$LM)
LMdata2_long$LM2 <- g2(LMdata2_long$LM)

LMdata2_long$group_int<-LMdata2_long$y*LMdata2_long$group
LMdata_long$group_int<-LMdata_long$y*LMdata_long$group

####################################
#Landmark supermodels
LMsupercox2 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 +
                       LM1 + LM2 + cluster(id), data=LMdata2_long, method="breslow")
bh_supercox2<-basehaz(LMsupercox2,centered=FALSE)

LMsupercox2_int <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ group + y + y_tau + y_tau2 +
                           group_int + 
                           LM1 + LM2 + cluster(id), data=LMdata2_long, method="breslow")
bh_supercox2_int<-basehaz(LMsupercox2_int,centered=FALSE)

###############################################################################
###############################################################################
### Copula Models
###############################################################################
###############################################################################
#Step 1: Model for T 
################################################
#Model for T 
cox_T<-coxph(Surv(Time,event)~group,data=train_data.id)

#Parametric Weibull hazard 
weib_T<-survreg(Surv(Time,event)~group,data=train_data.id,dist="weibull")
#   survreg's scale  =    1/(rweibull shape)
#   survreg's intercept = log(rweibull scale)
shape.weib<-1/weib_T$scale
scale.weib<-exp(weib_T$coef[1])
coef.weib<- -weib_T$coef[-1]/weib_T$scale

F_T_func_cond<-function(dat)
{
  xdata1<-as.numeric(dat[1])
  tt<-as.numeric(dat[2])
  tt_LM<-as.numeric(dat[3])
  
  bh<-list(basehaz(cox_T,FALSE))
  
  bet<-list(cox_T$coef)
  
  names(bet)<-c("cox_T")
  
  sfi<-bh
  sfi<-lapply(sfi,function(x) data.frame(x,haz0=diff(c(0,x$hazard))))
  
  sfi[[1]]["haz"]<-sfi[[1]]$haz0*as.numeric(exp(bet[["cox_T"]][1]*xdata1))
  
  sfi<-lapply(sfi,function(x) data.frame(x,Haz=cumsum(x$haz)))
  tmp<-lapply(sfi,function(x) evalstep(x$time,x$Haz,c(tt,tt_LM),subst=0))
  Fw<-lapply(tmp,function(x) ((1-exp(-x[1]))-(1-exp(-x[2])))/exp(-x[2]))
  
  Fw_weib_tt<-exp(-(tt/scale.weib)^shape.weib*exp(coef.weib*xdata1))
  Fw_weib_tt_LM<-exp(-(tt_LM/scale.weib)^shape.weib*exp(coef.weib*xdata1))
  Fw_weib<-((1-Fw_weib_tt)-(1-Fw_weib_tt_LM))/Fw_weib_tt_LM
  
  return(c(unlist(Fw),Fw_weib)) 
}

FT_cond_vec<-apply(LMdata_long[,c("group","Time","LM")],1,F_T_func_cond)
FT_cond_vec<-t(FT_cond_vec)
FT_names<-c("cox","weib")
FT_cond_vec<-data.frame(FT_cond_vec)
colnames(FT_cond_vec)<-FT_names
#FT_cond_vec contains the conditional survival predictions for patients in the training data set for all 
#of the survival models considered (different columns)

################################################
#Step 2: Model for Z (no longer need model for Z*!) 
################################################
F_Z<-function(z,mean_Z,sigma_Z,log.prob=FALSE)
{
  pnorm(z,mean=mean_Z,sd=sigma_Z,log.p=log.prob)
}

f_Z<-function(z,mean_Z,sigma_Z)
{
  dnorm(z,mean=mean_Z,sd=sigma_Z)
}

#Models for mean, sd of Z 
simp.red<-~LM+group-1
bspint.red<-~bs(LM, Boundary.knots=c(0,cens_horiz))*group-1

mods_Z<-list(simp.red,bspint.red) 
mods_names<-c("simp","bspint") 
names(mods_Z)<-mods_names

mean_names<-c("bspint")
sd_names<-c("simp")
             
#Fitting model for Z allowing mead and variance to be the functions considered above
#NOTE: choose which models you want to fit to minimize running time 
for(mean_string in mean_names)
{
  for(sd_string in sd_names)
  {
    mod_mean_Z<-mods_Z[[mean_string]]
    mod_sd_Z<-mods_Z[[sd_string]]
    
    xmat.mean<-model.matrix(mod_mean_Z,data=LMdata_long)
    xmat.var<-model.matrix(mod_sd_Z,data=LMdata_long)
    
    fit<-lmvar(LMdata_long$y,xmat.mean,xmat.var)
    assign(paste0("mod_Z_full_mean",mean_string,"_sd",sd_string),fit)
    assign(paste0("mod_Z_mean",mean_string,"_sd",sd_string),coef(fit))
  }
}


##################################
#Joint distribution
##################################
FZT<-function(FT_cond,z,mean_Z,sigma_Z,rho,cop,df)
{
  if(cop=="gaussian") {
    a<-qnorm(FT_cond,mean=0,sd=1)
    b<-qnorm(F_Z(z,mean_Z,sigma_Z),mean=0,sd=1)
    ret<-ifelse(a==-Inf,0,
                ifelse(a==Inf,F_Z(0,mean_Z,sigma_Z),
                       ifelse(b==Inf,FT_cond,pbivnorm(x=a,y=b,rho=rho))))
  } else if(cop=="t") {
    a<-qt(FT_cond,df=df)
    b<-qt(F_Zstar(zstar,mean_Zstar,sigma_Zstar),df=df)
    ret<-ifelse(a==-Inf,0,
                ifelse(a==Inf,F_Zstar(0,mean_Zstar,sigma_Zstar),
                       ifelse(b==Inf,FT_cond,pmvt(lower=c(-Inf,-Inf),upper=c(a,b),corr=rho*matrix(c(1,0,0,1),nrow=2,df=df)))))
  }
  return(ret)
}

##################################
#Maximize joint likelihood
##################################
#Association models
loglik_ZT<-function(g,FTcond,meanZ,sigmaZ,dat,y,copula="gaussian",v=NA)
{
  #L1: P(T=t,Z=z)=f(t,z)
  #L2: P(T>t,Z=z)
  
  eta<-coef_dat%*%g
  rho<-1-2/(1+exp(2*eta))
  
  if(copula=="gaussian") {
    q1<-qnorm(FTcond,mean=0,sd=1)
    q2<-qnorm(F_Z(y,meanZ,sigmaZ,log.prob=TRUE),mean=0,sd=1,log.p=TRUE)
    
    L1<--1/2*log(1-rho^2)-(rho^2*(q1^2+q2^2)-2*rho*q1*q2)/(2*(1-rho^2))
    L2<-log(pnorm(-(q1-rho*q2)/sqrt(1-rho^2)))
    
  } else if(copula=="t") {
    q1<-qt(FTcond,df=v)
    q2<-qt(F_Z(y,meanZ,sigmaZ,log.prob=TRUE),df=v,log.p=TRUE)
    
    L1<--(v+1)/2*log(1+1/(v*(1-rho^2))*(q1^2+q2^2-2*rho*q1*q2))-dt(q1,df=v,log=TRUE)-dt(q2,df=v,log=TRUE)-1/2*log(1-rho^2)
    L2<-log(pt(-(q1-rho*q2)/sqrt((v+q2^2)*(1-rho^2)/(v+1)),df=v+1))
  }
  
  L1_cont<-L1[which(dat$event==1)]
  L2_cont<-L2[which(dat$event==0)]
  
  ret<--(sum(L1_cont)+sum(L2_cont))
  
  return(ret)
}

#Perform maximization for the different copula models
#Error for one-dimensional optimization if "const" model considered 
simp<-~LM+group
bsp<-~bs(LM, Boundary.knots=c(0,cens_horiz))+group
int<-~LM*group

mods_rho<-list(simp,bsp,int) 
mods_rho_names<-c("simp","bsp","int") 
names(mods_rho)<-mods_rho_names

rho_names<-c("simp","bsp","int")

for(cop_string in c("gaussian","t")) 
{
  for(rho_string in rho_names)
  {
    for(mean_string in mean_names)
    {
      for(sd_string in sd_names)
      {
        for(FT_string in FT_names)
        {
          xmat.mean.pred<-model.matrix(mods_Z[[mean_string]],LMdata_long)
          xmat.sd.pred<-model.matrix(mods_Z[[sd_string]],LMdata_long)
          modZ_vals<-predict(get(paste0("mod_Z_full_mean",mean_string,"_sd",sd_string)),xmat.mean.pred,xmat.sd.pred)
          meanZ<-modZ_vals[,1]
          sdZ<-modZ_vals[,2]
          
          coef_dat<-model.matrix(mods_rho[[rho_string]],data=LMdata_long)
          
          start_val<-rep(0,ncol(coef_dat))
          temp<-try(optim(start_val,loglik_ZT,FTcond=FT_cond_vec[,FT_string],
                          meanZ=meanZ,
                          sigmaZ=sdZ,
                          y=LMdata_long$y,
                          v=df0,
                          copula=cop_string,
                          dat=LMdata_long,
                          method="Nelder-Mead",control=list(maxit=5000,trace=FALSE)))
          assign(paste0("MLE_rho",rho_string,"_mean",mean_string,"_sd",sd_string,"_",FT_string,"_",cop_string),temp$par) 
        }
      }
    }
  }
}
