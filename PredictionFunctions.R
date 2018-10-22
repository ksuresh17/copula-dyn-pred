#Prediction functions for Continuous biomarker

#######################################################
#Landmark Models 
#######################################################
#LM0 
BSpredict_LM0<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  Fw <- NULL
  tti<-t0
  for(i in 1:length(zdata))
  {
    sfi<-subset(bh,bh$strata==paste0("LM=",tti))
    sfi$haz0<-diff(c(0,sfi$hazard))
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i]+bet[2]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LMInt0
BSpredict_LMInt0<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  Fw <- NULL
  tti<-t0
  for(i in 1:length(zdata))
  {
    sfi<-subset(bh,bh$strata==paste0("LM=",tti))
    sfi$haz0<-diff(c(0,sfi$hazard))
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i]+bet[2]*zdata[i]+bet[3]*xdata[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM1 
BSpredict_LM1<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  nt <- length(zdata)
  tti<-t0
  Fw <- NULL
  for(i in 1:nt)
  {
    sfi<-subset(bh,bh$strata==paste0("LM=",tti))
    sfi$haz0<-diff(c(0,sfi$hazard))
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*(tti) + bet[4]*zdata[i]*(tti)^2)
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LMInt1
BSpredict_LMInt1<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  nt <- length(zdata)
  tti<-t0
  Fw <- NULL
  for(i in 1:nt)
  {
    sfi<-subset(bh,bh$strata==paste0("LM=",tti))
    sfi$haz0<-diff(c(0,sfi$hazard))
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*(tti) + bet[4]*zdata[i]*(tti)^2 + 
                                bet[5]*xdata[i]*zdata[i])
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LMInt2 
BSpredict_LM2<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  Fw <- NULL
  tti<-t0
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  for(i in 1:length(zdata))
  {
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*(tti) + bet[4]*zdata[i]*(tti) +
                                bet[5]*g1(tti) + bet[6]*g2(tti)) 
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}


#LMInt2
BSpredict_LMInt2<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  nt <- length(zdata)
  Fw <- NULL
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  for(i in 1:nt)
  {
    tti<-t0
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*(tti) + bet[4]*zdata[i]*(tti) +
                                bet[5]*xdata[i]*zdata[i] + 
                                bet[6]*g1(tti) + bet[7]*g2(tti)) 
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM3
BSpredict_LM3<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  nt <- length(zdata)
  Fw <- NULL
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  for(i in 1:nt)
  {
    tti<-t0
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*(sfi$time-tti) + bet[4]*zdata[i]*(sfi$time-tti)^2 +
                                bet[5]*g1(tti) + bet[6]*g2(tti))
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LMInt3
BSpredict_LMInt3<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  nt <- length(zdata)
  Fw <- NULL
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  for(i in 1:nt)
  {
    tti<-t0
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*(sfi$time-tti) + bet[4]*zdata[i]*(sfi$time-tti)^2 +
                                bet[5]*xdata[i]*zdata[i] + 
                                bet[6]*g1(tti) + bet[7]*g2(tti)) 
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LM4
BSpredict_LM4<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  nt <- length(zdata)
  Fw <- NULL
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  for(i in 1:nt)
  {
    tti<-t0
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*tti + bet[4]*zdata[i]*tti^2 + 
                                bet[5]*zdata[i]*(sfi$time-tti) + bet[6]*zdata[i]*(sfi$time-tti)^2 +
                                bet[7]*g1(tti) + bet[8]*g2(tti))
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#LMInt4
BSpredict_LMInt4<-function(bh,bet,zdata,t0,xdata,w_predict)
{
  nt <- length(zdata)
  Fw <- NULL
  sfi<-bh
  sfi$haz0<-diff(c(0,sfi$hazard))
  for(i in 1:nt)
  {
    tti<-t0
    sfi$hazard<-sfi$haz0* exp(bet[1]*xdata[i] + 
                                bet[2]*zdata[i] + bet[3]*zdata[i]*tti + bet[4]*zdata[i]*tti^2 + 
                                bet[5]*zdata[i]*(sfi$time-tti) + bet[6]*zdata[i]*(sfi$time-tti)^2 +
                                bet[7]*zdata[i]*xdata[i] + 
                                bet[8]*g1(tti) + bet[9]*g2(tti))
    sfi$Haz <- cumsum(sfi$hazard)
    tmp <- evalstep(sfi$time,sfi$Haz,c(tti,tti+w_predict),subst=0)
    Fw <- c(Fw,1-exp(-(tmp[2]-tmp[1])))
  }
  return(Fw)
}

#######################################################
#Joint Models 
#######################################################
BSpredict_JM<-function(JM_mod,sub_data,t0,w_predict,sim)
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

#######################################################
#True Model
#######################################################
BSpredict_Truth<-function(trueVals,sub_data,t0,w_predict)
{
  W <- cbind("(Intercept)" = 1, "Group" = sub_data$group)
  eta.t <- as.vector(W %*% trueVals$gammas)
  
  Fw<-NULL
  for(i in 1:nrow(sub_data))
  {
    h <- function (s) {
      group <- sub_data$group[i]
      XX <- cbind(1, s, group, group*s)
      ZZ <- cbind(1, s)
      f1 <- as.vector(XX %*% trueVals$betas + rowSums(ZZ * trueVals$b_test[rep(sub_data$id[i], nrow(ZZ)), ]))
      exp(log(trueVals$phi) + (trueVals$phi - 1) * log(s) + eta.t[i] + f1 * trueVals$alpha)
    }
    Fw<-c(Fw,1-exp(-integrate(h, lower = t0, upper = t0+w_predict)$value))
  }
  return(Fw)
}

########################################################
#Copula Prediction function 
########################################################
BSpredict_Copula<-function(zdata,t0,xdata,w_predict,mean_string,sd_string,F_T_dat,FT_string,rho_string,obs.time=t0) #,v)
{
  covdat_Z<-data.frame(LM=obs.time,group=xdata)
  
  covdat_rho<-data.frame(LM=t0,group=xdata)
  
  #modZ
  xmat.mean.pred<-model.matrix(mods.red[[mean_string]],covdat_Z)
  xmat.sd.pred<-model.matrix(mods.red[[sd_string]],covdat_Z)
  modZ_vals<-predict(get(paste0("mod_Zstar_full_mean",mean_string,"_sd",sd_string)),xmat.mean.pred,xmat.sd.pred)
  meanZ<-modZ_vals[,1]
  sdZ<-modZ_vals[,2]
  
  rho_param<-get(paste0("MLE_rho",rho_string,"_mean",mean_string,"_sd",sd_string,"_",FT_string))
  
  coef_dat<-model.matrix(mods[[rho_string]],data=covdat_rho)
  
  eta_Ystar_LM<-coef_dat%*%rho_param
  rho_LM<-1-2/(1+exp(2*eta_Ystar_LM))
  
  fZ_pred<-f_Zstar(zstar=zdata,mean=meanZ,sigma_Zstar=sdZ)
  FZ_pred<-F_Zstar(zstar=zdata,mean=meanZ,sigma_Zstar=sdZ,log.prob=TRUE)
  
  FT_tau_w_cond<-F_T_dat[,FT_string]
  
  q2_pred<-qnorm(FZ_pred,log.p=TRUE)
  
  pred<-apply(FT_tau_w_cond,2,function(x) {
    q1_pred<-qnorm(x)
    norm_pred<-pnorm((q1_pred-rho_LM*q2_pred)/sqrt(1-rho_LM^2))
    return(norm_pred)
  })
  
  return(pred)
}

#######################################################
#Code obtianed from Blanche et al. (2015) for computation of dynamic BS and AUC
#######################################################
#Brier Score function from Blanche et al. (2015)
#PE function from Blanche et al. (2015)
PE<-function(model,LMx,w_predict,BS_dat)
{    
  AUCst.M1 <- rep(NA,length(LMx))
  BrierS.s.M1 <- rep(NA,length(LMx))
  for (s in LMx){
    # print(s)
    # Create landmark data set
    d.s<-BS_dat[,c("Time","event")]
    d.s$Pred.s.M1<-BS_dat[,model]
    d.s<-subset(d.s,BS_dat$LM==s) #subset(d.s,d.s$Time>s)
    d.s$time.s<-d.s$Time-s
    # AUC and BS for prediction based on M1
    # estimate ROC curve and AUC
    ROC.s.M1<-timeROC(T=d.s$time.s,
                      delta=d.s$event,
                      marker=d.s$Pred.s.M1,
                      cause=1,weighting="marginal",
                      times=c(w_predict),
                      iid=TRUE)
    # estimate expected Brier score
    BS.s.M1 <- BS(timepoints=c(w_predict),
                  times=d.s$time.s,
                  status=d.s$event,
                  pred=as.matrix(d.s$Pred.s.M1),
                  cause=1)
    # save useful results (estimates, s.e. and iid decompositions)
    BrierS.s.M1[which(s==LMx)] <- BS.s.M1$BS # BS estimate
    AUCst.M1[which(s==LMx)] <- ROC.s.M1$AUC[2] # AUC estimate
  }
  return(list(BrierS.s.M1,AUCst.M1))
}

BS <- function(timepoints,times,status,pred,cause=1){ 
  n <- length(times)
  n_times <- length(timepoints)
  timepoints <- timepoints[order(timepoints)]
  times_names <- paste("t=",timepoints,sep="")
  # output initialisation 
  BS <- rep(NA,n_times)
  CumInci <- rep(NA,n_times)
  surv <- rep(NA,n_times)
  Stats <- matrix(NA,nrow=n_times,ncol=4)
  hit1_all <- matrix(NA,nrow=n,ncol=n_times)
  hit2_all <- matrix(NA,nrow=n,ncol=n_times)
  epsilon_i <- matrix(NA,nrow=n,ncol=n_times)
  #adds name to outputs
  names(BS) <- times_names
  names(CumInci) <- times_names
  names(surv) <- times_names
  colnames(Stats) <- c("Cases","survivor at t","Other events at t","Censored at t")
  rownames(Stats) <- times_names
  colnames(epsilon_i) <- times_names
  colnames(hit1_all) <-  times_names
  colnames(hit2_all)  <- times_names 
  # we need to order to use the ipcw() function of the pec package
  #browser()
  order_T <- order(times)
  times <-  times[order_T]
  delta  <-  status[order_T]
  pred <-  pred[order_T,,drop=FALSE]
  #compute KM weights
  weights <- ipcw(Surv(failure_time,status)~1,
                  data=data.frame(failure_time=times,status=as.numeric(delta!=0)),
                  method="marginal",times=timepoints,subjectTimes=times,subjectTimesLag=1)
  Mat_data <- cbind(times,delta,as.numeric(delta==0))
  colnames(Mat_data) <- c("T","delta","indic_Cens")
  # computate weights of cases
  Weights_cases_all <- 1/(weights$IPCW.subjectTimes*n)
  # loop on all time points
  for(t in 1:n_times){
    Cases <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]==cause)
    Controls_1 <- (Mat_data[,"T"]> timepoints[t] )
    Controls_2 <- (Mat_data[,"T"]<= timepoints[t] &  Mat_data[,"delta"]!=cause & Mat_data[,"delta"]!=0)  
    # compute weights
    Weights_controls_1 <- rep(1/(weights$IPCW.times[t]*n),times=n)
    Weights_cases <- Weights_cases_all
    Weights_controls_2 <- Weights_cases_all
    Weights_cases[!Cases] <- 0
    Weights_controls_1[!Controls_1] <- 0
    Weights_controls_2[!Controls_2] <- 0   
    #compute outputs
    CumInci[t] <- c(sum(Weights_cases))
    surv[t] <- c(sum(Weights_controls_1))
    Stats[t,] <- c(sum(Cases),sum(Controls_1),sum(Controls_2),n-sum(Cases)-sum(Controls_1)-sum(Controls_2)) 
    hit1_all[,t] <- (Weights_controls_1*((pred[,t])^2))*n
    hit2_all[,t] <- (Weights_cases*((1-pred[,t])^2) + Weights_controls_2*((pred[,t])^2))*n
    BS[t] <- (sum(hit1_all[,t]) +sum(hit2_all[,t]))/n
  } 
  
  out <- list(BS=BS,res=(hit1_all+hit2_all),
              CumulativeIncidence=CumInci,survProb=surv,n=n,Stats=Stats,timepoints=timepoints
  )
  class(out) <- "ipcwEBS"
  out 
}
