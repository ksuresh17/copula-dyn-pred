########################
#Code for application of Copula method to heart valve data 
########################
rm(list=ls())
library(SmoothHazard)
library(dynpred)
library(lava)
library(prodlim)
library(mstate)
library(reshape2)
library(ggplot2)
library(cowplot)
library(plyr)
library(pbivnorm)
library(splines)
library(lmvar)

library(JM)
library(joineR)

data(heart.valve)

heart.valve<-subset(heart.valve,!is.na(heart.valve$log.lvmi))
heart.valve<-subset(heart.valve,fuyrs>time)
heart.valve$hs<-as.numeric(heart.valve$hs)-1
#set prediction window
w<-4 
cens_horiz<-max(heart.valve$time)

#Order heartvalve data
heart.valve<-heart.valve[order(heart.valve$num,heart.valve$time),]
llvmi0<-heart.valve[!duplicated(heart.valve$num),]
heart.valve<-merge(x=heart.valve,y=llvmi0[,c("num","log.lvmi")],by="num",suffixes=c("","0"))

#Create long dataset for use in fitting landmark models 
dat.id<-unique(heart.valve[,c("fuyrs","status","age","hs","sex","num","size","bsa","log.lvmi0")])
long_data<-heart.valve
long_data$LM<-long_data$time
long_data$survtime<-long_data$fuyrs
long_data$wsurvtime<-pmin(long_data$LM+w,long_data$survtime)
long_data$death<-long_data$status
long_data$survstatus<-ifelse(long_data$survtime<=long_data$wsurvtime&long_data$death==1,1,0)
long_data$inspection.time<-long_data$LM
long_data<-subset(long_data,long_data$LM<long_data$wsurvtime)
long_data$illness<-long_data$log.lvmi
long_data$id<-long_data$num

long_data<-subset(long_data,long_data$LM<=8)
LMdata<-long_data
LMdata<-LMdata[order(LMdata$id,LMdata$inspection.time),]
tt<-sort(unique(LMdata$wsurvtime[LMdata$survstatus==1]))

LMdata$Tstart<-LMdata$inspection.time
LMdata2<-survSplit(Surv(Tstart,wsurvtime,survstatus)~.,data=LMdata,cut=tt,end="wsurvtime",start="Tstart",event="survstatus")

#Define covariates used in landmark models 
LMdata2$illnessX<-LMdata2$illness*(LMdata2$wsurvtime-LMdata2$inspection.time)
LMdata2$illnessX2<-LMdata2$illness*(LMdata2$wsurvtime-LMdata2$inspection.time)^2
LMdata2$illness_tau<-LMdata2$illness*(LMdata2$inspection.time)
LMdata2$illness_tau2<-LMdata2$illness*(LMdata2$inspection.time)^2

g1<-function(t) t
g2<-function(t) t^2
LMdata2$LM1<-g1(LMdata2$LM)
LMdata2$LM2<-g2(LMdata2$LM)

LMdata2<-subset(LMdata2,LMdata2$wsurvtime-LMdata2$Tstart>.0001)

####################################################################
#Landmark Models
####################################################################
#LM2 
LMsupercox2 <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ illness + illness_tau + illness_tau2 + LM1+LM2+ 
                       age + hs +sex +cluster(id), data=LMdata2, method="breslow")

bh_supercox2<-basehaz(LMsupercox2,centered=FALSE)

#LM3
LMsupercox3<-coxph(Surv(Tstart,wsurvtime,survstatus)~illness+illnessX+illnessX2+LM1+LM2+
                     age + hs +sex +cluster(id),data=LMdata2,method="breslow")

bh_supercox3<-basehaz(LMsupercox3,centered=FALSE)

#LM2Int 
LMsupercox2Int <- coxph(Surv(Tstart,wsurvtime,survstatus) ~ illness + illness_tau + illness_tau2 + LM1+LM2+ 
                          age + hs +sex + 
                          age*illness + hs*illness + sex*illness+ 
                          cluster(id), data=LMdata2, method="breslow")

bh_supercox2Int<-basehaz(LMsupercox2Int,centered=FALSE)

#LM3Int
LMsupercox3Int<-coxph(Surv(Tstart,wsurvtime,survstatus)~illness+illnessX+illnessX2+LM1+LM2+
                        age + hs +sex +
                        age*illness + hs*illness + sex*illness+
                        cluster(id),data=LMdata2,method="breslow")

bh_supercox3Int<-basehaz(LMsupercox3Int,centered=FALSE)

#Prediction function for LMInt3
BSpredict_LMInt3<-function(bh,bet,zdata,t0,age_data,sex_data,hs_data,w_predict)
{
  nt <- length(zdata)
  Fw <- NULL
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  for(i in 1:nt)
  {
    tti<-t0
    sfi$hazard<-sfi$haz0* exp(bet[1]*zdata[i] + bet[2]*zdata[i]*(sfi$time-tti) + bet[3]*zdata[i]*(sfi$time-tti)^2 +
                                bet[4]*g1(tti) + bet[5]*g2(tti)+
                                bet[6]*age_data[i]+bet[7]*hs_data[i]+bet[8]*sex_data[i]+
                                bet[9]*age_data[i]*zdata[i]+bet[10]*hs_data[i]*zdata[i]+sex_data[i]*zdata[i]) 
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

####################################################################
#Joint Models
####################################################################
lmeFit<-lme(log.lvmi~time+age+hs+sex,data=heart.valve,random=~1|num)
survFit <- coxph(Surv(fuyrs, status) ~ age+hs+sex, data = dat.id,
                 x = TRUE)
jointFit1 <- jointModel(lmeFit, survFit, timeVar = "time") #,method="Cox-PH-GH")

lmeFit2<-lme(log.lvmi~time+age+hs+sex,data=heart.valve,random=~time|num)
survFit2 <- coxph(Surv(fuyrs, status) ~ age+hs+sex, data = dat.id,
                  x = TRUE)
jointFit2 <- jointModel(lmeFit2, survFit2, timeVar = "time") #,method="Cox-PH-GH")

lmeFit3<-lme(log.lvmi~time+age+hs+sex+time:hs+time:sex,data=heart.valve,random=~time|num)
survFit3 <- coxph(Surv(fuyrs, status) ~ age+hs+sex, data = dat.id,
                  x = TRUE)
jointFit3 <- jointModel(lmeFit3, survFit3, timeVar = "time") #,method="Cox-PH-GH")

#Prediction function for JMs
BSpredict_JM<-function(JM_mod,sub_data,t0,w_predict,sim=TRUE)
{
  pred_cond<-survfitJM(JM_mod,
                       newdata=sub_data,
                       idVar="id",survTimes = t0+w_predict,simulate=sim)
  if(sim==FALSE)
  {
    Fw<-1-t(matrix(unlist(pred_cond$summaries),nrow=2))[,2]
  } else 
  {
    Fw<-1-t(matrix(unlist(pred_cond$summaries),nrow=5))[,2]
  }
  return(Fw)
}

##################################################################################################
#COPULA
##################################################################################################
#Step 1: Model for T 
#Fit cox_T, weib_T
################################################
#Model for T 
cox_T<-coxph(Surv(fuyrs,status)~age+hs+sex,data=dat.id)
#Test nonPH
zp<-cox.zph(cox_T)

weib_T<-survreg(Surv(fuyrs,status)~age+hs+sex,data=dat.id,dist="weibull")
#   survreg's scale  =    1/(rweibull shape)
#   survreg's intercept = log(rweibull scale)
shape.weib<-1/weib_T$scale
scale.weib<-exp(weib_T$coef[1])
coef.weib<- -weib_T$coef[-1]/weib_T$scale

#Function that returns F values for all cox models considered
F_T_func_cond<-function(dat)
{
  xdata1<-as.numeric(dat[1:3])
  tt<-as.numeric(dat[4])
  tt_LM<-as.numeric(dat[5])
  
  bh<-list(basehaz(cox_T,FALSE))
  
  bet<-list(cox_T$coef)
  
  names(bet)<-c("cox_T")
  
  sfi<-bh
  sfi<-lapply(sfi,function(x) data.frame(x,haz0=diff(c(0,x$hazard))))
  
  sfi[[1]]["haz"]<-sfi[[1]]$haz0*as.numeric(exp(bet[["cox_T"]]%*%xdata1))
  
  sfi<-lapply(sfi,function(x) data.frame(x,Haz=cumsum(x$haz)))
  tmp<-lapply(sfi,function(x) evalstep(x$time,x$Haz,c(tt,tt_LM),subst=0))
  Fw<-lapply(tmp,function(x) ((1-exp(-x[1]))-(1-exp(-x[2])))/exp(-x[2]))  #-exp(-x[1])/exp(-x[2])+1) #(-exp(-x[1])+exp(-x[2]))/(exp(-x[2])))
  
  Fw_weib_tt<-exp(-(tt/scale.weib)^shape.weib*exp(coef.weib%*%xdata1))
  Fw_weib_tt_LM<-exp(-(tt_LM/scale.weib)^shape.weib*exp(coef.weib%*%xdata1))
  Fw_weib<-((1-Fw_weib_tt)-(1-Fw_weib_tt_LM))/Fw_weib_tt_LM
  
  return(c(unlist(Fw),Fw_weib))
}

FT_cond_vec<-apply(LMdata[,c("age","hs","sex","fuyrs","LM")],1,F_T_func_cond)
FT_cond_vec<-t(FT_cond_vec)
FT_names<-c("cox", "weib") 
colnames(FT_cond_vec)<-FT_names

################################################
#Step 2: Model for Z
################################################
F_Zstar<-function(zstar,mean_Zstar,sigma_Zstar,log.prob=FALSE)
{
  pnorm(zstar,mean=mean_Zstar,sd=sigma_Zstar,log.p=log.prob)
}

f_Zstar<-function(zstar,mean_Zstar,sigma_Zstar)
{
  dnorm(zstar,mean=mean_Zstar,sd=sigma_Zstar)
}

#Models for Z 
const<-~1
fLM<-~LM
covs<-~age+hs+sex
simp<-~LM+age+hs+sex
int<-~LM*(age+hs+sex)
sq<-~LM+I(LM^2)+age+hs+sex
sqint<-~LM*(age+hs+sex)+I(LM^2)*(age+hs+sex)
sqpartint<-~LM*age+LM*hs+LM*sex+I(LM^2)
bsp<-~bs(LM, Boundary.knots=c(0,cens_horiz))+age+hs+sex
bspint<-~bs(LM, Boundary.knots=c(0,cens_horiz))*(age+hs+sex)

mods<-list(const,fLM,covs,simp,int,sq,sqint,sqpartint,bsp,bspint)

mods_names<-c("const","fLM","covs","simp","int","sq","sqint","sqpartint","bsp","bspint")
names(mods)<-mods_names

#Reduced models for Z 
const.red<-~1-1
fLM.red<-~LM-1
covs.red<-~age+hs+sex-1
simp.red<-~LM+age+hs+sex-1
int.red<-~LM*(age+hs+sex)-1
sq.red<-~LM+I(LM^2)+age+hs+sex-1
sqint.red<-~LM*age+I(LM^2)*age+LM*hs+I(LM^2)*hs+LM*sex+I(LM^2)*sex-1
sqpartint.red<-~LM*age+LM*hs+LM*sex+I(LM^2)-1
bsp.red<-~bs(LM, Boundary.knots=c(0,cens_horiz))+age+hs+sex-1
bspint.red<-~bs(LM, Boundary.knots=c(0,cens_horiz))*(age+hs+sex)-1

mods.red<-list(const.red,fLM.red,covs.red,simp.red,int.red,sq.red,sqint.red,sqpartint.red,bsp.red,bspint.red)
names(mods.red)<-mods_names

for(mean_string in mods_names)
{
  for(sd_string in mods_names)
  {
    print(c(mean_string,sd_string))
    mod_mean_Zstar<-mods.red[[mean_string]]
    mod_sd_Zstar<-mods.red[[sd_string]]
    
    xmat.mean<-model.matrix(mod_mean_Zstar,data=LMdata)
    xmat.var<-model.matrix(mod_sd_Zstar,data=LMdata)
    
    fit<-lmvar(LMdata$log.lvmi,xmat.mean,xmat.var)
    assign(paste0("mod_Zstar_full_mean",mean_string,"_sd",sd_string),fit)
    assign(paste0("mod_Zstar_mean",mean_string,"_sd",sd_string),coef(fit))
  }
}

##################################
#Joint distribution
##################################
FZstarT<-function(FT_cond,zstar,mean_Zstar,sigma_Zstar,rho,df)
{
  a<-qnorm(FT_cond,mean=0,sd=1)
  b<-qnorm(F_Zstar(zstar,mean_Zstar,sigma_Zstar),mean=0,sd=1)
  ret<-ifelse(a==-Inf,0,
              ifelse(a==Inf,F_Zstar(0,mean_Zstar,sigma_Zstar),
                     ifelse(b==Inf,FT_cond,pbivnorm(x=a,y=b,rho=rho))))
  return(ret)
}

##################################
#Maximize joint likelihood
##################################
#Association models
loglik_ZT<-function(g,FTcond,meanZ,sigmaZ,dat) #,v)
{
  eta<-coef_dat%*%g
  rho<-1-2/(1+exp(2*eta))
  
  q1<-qnorm(FTcond,mean=0,sd=1)
  q2<-qnorm(F_Zstar(dat$illness,meanZ,sigmaZ,log.prob=TRUE),mean=0,sd=1,log.p=TRUE)
  
  #L1: P(T=t,Z=z)=f(t,z)
  L1<--1/2*log(1-rho^2)-(rho^2*(q1^2+q2^2)-2*rho*q1*q2)/(2*(1-rho^2))
  L1_cont<-L1[which(dat$status==1)]
  
  #L2: P(T>t,Z=z)
  L2<-log(pnorm(-(q1-rho*q2)/sqrt(1-rho^2)))
  L2_cont<-L2[which(dat$status==0)]
  
  ret<--(sum(L1_cont)+sum(L2_cont))
  
  return(ret)
}

rho_names<-c("simp")
mean_names<-c("bsp")
sd_names<-c("const")

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
                        dat=LMdata,
                        method="BFGS",control=list(maxit=5000,trace=FALSE))) #,hessian=TRUE)
        
        assign(paste0("MLE_rho",rho_string,"_mean",mean_string,"_sd",sd_string,"_",FT_string),temp$par)
      }
    }
  }
}
########################################################
#Prediction function for Copula
########################################################
BSpredict_Copula<-function(zdata,t0,age_data,hs_data,sex_data,w_predict,mean_string,sd_string,FT_string,rho_string,obs.time=t0) 
{
  Fw<-NULL
  
  for(i in 1:length(zdata))
  {
    covdat_Z<-data.frame(LM=obs.time[i],age=age_data[i],hs=hs_data[i],sex=sex_data[i])
    
    covdat_rho<-data.frame(LM=t0,age=age_data[i],hs=hs_data[i],sex=sex_data[i])
    
    #modZ
    xmat.mean.pred<-model.matrix(mods.red[[mean_string]],rbind(covdat_Z,covdat_Z))
    xmat.sd.pred<-model.matrix(mods.red[[sd_string]],rbind(covdat_Z,covdat_Z))
    modZ_vals<-predict(get(paste0("mod_Zstar_full_mean",mean_string,"_sd",sd_string)),xmat.mean.pred,xmat.sd.pred)
    meanZ<-modZ_vals[1,1]
    sdZ<-modZ_vals[1,2]
    
    rho_param<-get(paste0("MLE_rho",rho_string,"_mean",mean_string,"_sd",sd_string,"_",FT_string))
    
    coef_dat<-model.matrix(mods[[rho_string]],data=covdat_rho)
    
    eta_Ystar_LM<-coef_dat%*%rho_param
    rho_LM<-1-2/(1+exp(2*eta_Ystar_LM))
    
    fZ_pred<-f_Zstar(zstar=zdata[i],mean=meanZ,sigma_Zstar=sdZ)
    FZ_pred<-F_Zstar(zstar=zdata[i],mean=meanZ,sigma_Zstar=sdZ,log.prob=TRUE)
    
    F_T_temp<-F_T_func_cond(c(age_data[i],hs_data[i],sex_data[i],t0+w_predict,t0))
    names(F_T_temp)<-FT_names
    FT_tau_w_cond<-F_T_temp[FT_string]
    
    q2_pred<-qnorm(FZ_pred,log.p=TRUE)
    q1_pred<-qnorm(FT_tau_w_cond)
    pred<-pnorm((q1_pred-rho_LM*q2_pred)/sqrt(1-rho_LM^2))
    
    Fw<-c(Fw,pred)
  }
  return(Fw)
}

BSpredict_Copula_Vec<-Vectorize(BSpredict_Copula,"w_predict")
BSpredict_JM_Vec<-Vectorize(BSpredict_JM,"w_predict")
BSpredict_LMInt3_Vec<-Vectorize(BSpredict_LMInt3,"w_predict")

#########################################################
#Figure creation (Comparison of methods)
#########################################################
#Select people still alive at LM1
LM_pred<-1
sub_LM1<-subset(LMdata,LMdata$survtime>LM_pred&LMdata$time<=LM_pred)
#Take last available marker measurement
sub2_LM1<-aggregate(sub_LM1, by = list(sub_LM1$id), FUN = tail, n = 1)
#Find quantiles of lvmi at that time
qlvmi<-quantile(sub2_LM1$log.lvmi)
# 0%      25%      50%      75%     100% 
# 4.244344 4.785697 4.999709 5.257337 6.351775 
sub2_LM1$cutlvmi<-cut(sub2_LM1$log.lvmi,c(0,5,10),include.lowest = TRUE)
qlvmi<-c(4.8,5.3)

palette_cols<-c("blue"='#377eb8',"red"='#e41a1c',"green"='#74c476')
par(mfrow=c(2,2))
#Plots of KM curves 
#Male, Stentless, llvmi<5
plot.times<-seq(0.1,10,.5)
i=1 
sex_ind=0
hs_ind=0
plot.dat<-subset(sub2_LM1,hs==hs_ind&cutlvmi==levels(sub2_LM1$cutlvmi)[i]&sub2_LM1$sex==sex_ind)
age_ind=median(plot.dat$age)
Cop_pred<-1-BSpredict_Copula_Vec(zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,hs_data=hs_ind,sex_data=sex_ind,
                                 w_predict=plot.times,mean_string="bsp",sd_string="const",FT_string="weib",rho_string="simp") 
JMsim_pred=1-BSpredict_JM_Vec(jointFit2,sub_data=data.frame(age=age_ind,sex=sex_ind,hs=hs_ind,id=301,log.lvmi=qlvmi[i],time=LM_pred),t0=LM_pred,w_predict=plot.times,sim=TRUE)
LMInt3_pred=1-BSpredict_LMInt3_Vec(bh=bh_supercox3Int,bet=coef(LMsupercox3Int),zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,sex_data=sex_ind,
                                   hs_data=hs_ind,w_predict=plot.times)

plot(survfit(Surv(fuyrs,status)~1,data=plot.dat),lty=1,conf.int=FALSE,xlim=c(0,10),lwd=2,xlab="Time (years)",ylab="Survival Probability",
     main="Male, Stentless, log(LVMI)<5",ylim=c(0,1))
lines(x=plot.times+LM_pred,y=Cop_pred,lty=2,col=palette_cols[1],lwd=2)
lines(x=plot.times+LM_pred,y=JMsim_pred,lty=3,col=palette_cols[2],lwd=2)
lines(x=plot.times+LM_pred,y=LMInt3_pred,lty=4,col=palette_cols[3],lwd=2)
abline(v=1,lty=2)
legend("bottomright",c("Kaplan-Meier","Copula","Joint Model","Landmark Model"),col=c("black",palette_cols[1:3]),lwd=2,lty=1:4)

#Male, Homograft, llvmi<5
i=1
sex_ind=0
hs_ind=1
plot.dat<-subset(sub2_LM1,hs==hs_ind&cutlvmi==levels(sub2_LM1$cutlvmi)[i]&sub2_LM1$sex==sex_ind)
age_ind=median(plot.dat$age)
Cop_pred<-1-BSpredict_Copula_Vec(zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,hs_data=hs_ind,sex_data=sex_ind,
                                 w_predict=plot.times,mean_string="bsp",sd_string="const",FT_string="weib",rho_string="simp") 
JMsim_pred=1-BSpredict_JM_Vec(jointFit2,sub_data=data.frame(age=age_ind,sex=sex_ind,hs=hs_ind,id=301,log.lvmi=qlvmi[i],time=LM_pred),t0=LM_pred,w_predict=plot.times,sim=TRUE)
LMInt3_pred=1-BSpredict_LMInt3_Vec(bh=bh_supercox3Int,bet=coef(LMsupercox3Int),zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,sex_data=sex_ind,
                                   hs_data=hs_ind,w_predict=plot.times)

plot(survfit(Surv(fuyrs,status)~1,data=plot.dat),lty=1,conf.int=FALSE,xlim=c(0,10),lwd=2,xlab="Time (years)",ylab="Survival Probability",
     main="Male, Homograft, log(LVMI)<5",ylim=c(0,1))
lines(x=plot.times+LM_pred,y=Cop_pred,lty=2,col=palette_cols[1],lwd=2)
lines(x=plot.times+LM_pred,y=JMsim_pred,lty=3,col=palette_cols[2],lwd=2)
lines(x=plot.times+LM_pred,y=LMInt3_pred,lty=4,col=palette_cols[3],lwd=2)
abline(v=1,lty=2)

#Male, Stentless, llvmi>5
i=2 
sex_ind=0
hs_ind=0
plot.dat<-subset(sub2_LM1,hs==hs_ind&cutlvmi==levels(sub2_LM1$cutlvmi)[i]&sub2_LM1$sex==sex_ind)
age_ind=median(plot.dat$age)
Cop_pred<-1-BSpredict_Copula_Vec(zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,hs_data=hs_ind,sex_data=sex_ind,
                                 w_predict=plot.times,mean_string="bsp",sd_string="const",FT_string="weib",rho_string="simp") 
JMsim_pred=1-BSpredict_JM_Vec(jointFit2,sub_data=data.frame(age=age_ind,sex=sex_ind,hs=hs_ind,id=301,log.lvmi=qlvmi[i],time=LM_pred),t0=LM_pred,w_predict=plot.times,sim=TRUE)
LMInt3_pred=1-BSpredict_LMInt3_Vec(bh=bh_supercox3Int,bet=coef(LMsupercox3Int),zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,sex_data=sex_ind,
                                   hs_data=hs_ind,w_predict=plot.times)

plot(survfit(Surv(fuyrs,status)~1,data=plot.dat),lty=1,conf.int=FALSE,xlim=c(0,10),lwd=2,xlab="Time (years)",ylab="Survival Probability",
     main="Male, Stentless, log(LVMI)>5")
lines(x=plot.times+LM_pred,y=Cop_pred,lty=2,col=palette_cols[1],lwd=2)
lines(x=plot.times+LM_pred,y=JMsim_pred,lty=3,col=palette_cols[2],lwd=2)
lines(x=plot.times+LM_pred,y=LMInt3_pred,lty=4,col=palette_cols[3],lwd=2)
abline(v=1,lty=2)

#Male, Homograft, llvmi>5
i=2
sex_ind=0
hs_ind=1
plot.dat<-subset(sub2_LM1,hs==hs_ind&cutlvmi==levels(sub2_LM1$cutlvmi)[i]&sub2_LM1$sex==sex_ind)
age_ind=median(plot.dat$age)
Cop_pred<-1-BSpredict_Copula_Vec(zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,hs_data=hs_ind,sex_data=sex_ind,
                                 w_predict=plot.times,mean_string="bsp",sd_string="const",FT_string="weib",rho_string="simp") 
JMsim_pred=1-BSpredict_JM_Vec(jointFit2,sub_data=data.frame(age=age_ind,sex=sex_ind,hs=hs_ind,id=301,log.lvmi=qlvmi[i],time=LM_pred),t0=LM_pred,w_predict=plot.times,sim=TRUE)
LMInt3_pred=1-BSpredict_LMInt3_Vec(bh=bh_supercox3Int,bet=coef(LMsupercox3Int),zdata=qlvmi[i],t0=LM_pred,age_data=age_ind,sex_data=sex_ind,
                                   hs_data=hs_ind,w_predict=plot.times)

plot(survfit(Surv(fuyrs,status)~1,data=plot.dat),lty=1,conf.int=FALSE,xlim=c(0,10),lwd=2,xlab="Time (years)",ylab="Survival Probability",
     main="Male, Homograft, log(LVMI)>5")
lines(x=plot.times+LM_pred,y=Cop_pred,lty=2,col=palette_cols[1],lwd=2)
lines(x=plot.times+LM_pred,y=JMsim_pred,lty=3,col=palette_cols[2],lwd=2)
lines(x=plot.times+LM_pred,y=LMInt3_pred,lty=4,col=palette_cols[3],lwd=2)
abline(v=1,lty=2)

#########################################################
#Figure creation: Predicted probabilities for two patients
#########################################################
#Individual A: Stentless, age 72, Male, no Death
p56<-subset(long_data,num==56)
#Individual B: Homograft, age 78, Male, Death
p76<-subset(long_data,num==76)

#Create plot of survival curve by time and log.lvmi observation 
#Individual A
par(mfrow=c(1,4))
sub.plot<-p56
for(i in 1:4)
{  
  plot.times<-seq(0,3,by=.1)
  par(mar=c(5,5,2,5))
  plot(x=sub.plot$time[1:i],y=sub.plot$log.lvmi[1:i],pch=16,ylim=c(4.8,5.2),xlab="Time (years)",ylab="log(LVMI)",xlim=c(0,8),cex.axis=1.2,cex=1.5,cex.lab=1.2)
  BS_pred<-1-BSpredict_Copula_Vec(zdata=sub.plot$log.lvmi[i],t0=sub.plot$time[i],age_data=sub.plot$age,hs_data=sub.plot$hs,sex_data=sub.plot$sex,
                                  w_predict=plot.times,mean_string="bsp",sd_string="const",FT_string="cox",rho_string="simp")
  
  par(new=TRUE)
  plot(x=c(0,plot.times)+sub.plot$time[i],c(1,BS_pred),type="l",ylim=c(0,1),xlim=c(0,8),xlab=NA,ylab=NA,axes=FALSE,lwd=2,col=cols[1])
  
  abline(v=sub.plot$time[i],lty=2)
  axis(side=4,las=2, cex.axis=1.2)
  mtext(side=4,line=3,"Survival Probability",cex=0.8)
}

#Individual B
par(mfrow=c(1,4))
sub.plot<-p76
for(i in 1:4)
{  
  plot.times<-seq(0,3,by=.1)
  par(mar=c(5,5,2,5))
  plot(x=sub.plot$time[1:i],y=sub.plot$log.lvmi[1:i],pch=16,ylim=c(5,5.5),xlab="Time (years)",ylab="log(LVMI)",xlim=c(0,8),cex.axis=1.2,cex=1.5,cex.lab=1.2)
  BS_pred<-1-BSpredict_Copula_Vec(zdata=sub.plot$log.lvmi[i],t0=sub.plot$time[i],age_data=sub.plot$age,hs_data=sub.plot$hs,sex_data=sub.plot$sex,
                                  w_predict=plot.times,mean_string="bsp",sd_string="const",FT_string="cox",rho_string="simp")
  
  par(new=TRUE)
  plot(x=c(0,plot.times)+sub.plot$time[i],c(1,BS_pred),type="l",ylim=c(0,1),xlim=c(0,8),xlab=NA,ylab=NA,axes=FALSE,lwd=2,col=palette_cols[1])
  
  abline(v=sub.plot$time[i],lty=2)
  abline(v=sub.plot$fuyrs[i],lwd=2,lty=3,col="blue")
  axis(side=4,las=2, cex.axis=1.2)
  mtext(side=4,line=3,"Survival Probability",cex=1)
}
