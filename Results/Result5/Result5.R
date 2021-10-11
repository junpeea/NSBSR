########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Script name : Result5

# Date : 201908xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################

rm(list=ls())
require(mcmc)
require(coda)
require(dplyr)
setwd("D:/NSBSR_FINAL")

Result5outMat = matrix(NA,nrow=12,10)
nd <- 22; ND <- 40

if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
  load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[1,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[1,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))

    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Result5outMat[1,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
    Result5outMat[1,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Result5outMat[1,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
    Result5outMat[1,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Result5outMat[1,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
    Result5outMat[1,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Result5outMat[1,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
    Result5outMat[1,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190926_N01m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N01m10.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[2,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[2,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists("THESIS_Results(New)/Result5/190926_P01m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P01m10.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[2,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[2,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[2,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[2,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[2,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[2,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[2,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[2,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190926_N01m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N01m25.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[3,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[3,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists("THESIS_Results(New)/Result5/190926_P01m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P01m25.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[3,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[3,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[3,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[3,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[3,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[3,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[3,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[3,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190926_N01m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N01m50.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[4,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[4,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists("THESIS_Results(New)/Result5/190926_P01m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P01m50.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[4,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[4,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[4,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[4,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[4,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[4,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[4,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[4,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}





if(file.exists("2019-08-05_NSBSR_Matern05_Grid(32,32)_ND=21_nd=21.Rdata")){
  load(file="2019-08-05_NSBSR_Matern05_Grid(32,32)_ND=21_nd=21.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
Result5outMat[5,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[5,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists("2019-08-05_PBSR_Matern05_Grid(32,32)_ND=21_nd=21.Rdata")){
  load(file="2019-08-05_PBSR_Matern05_Grid(32,32)_ND=21_nd=21.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[5,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[5,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[5,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[5,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[5,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[5,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[5,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[5,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190926_N05m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N05m10.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[6,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[6,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists("THESIS_Results(New)/Result5/190926_P05m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P05m10.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[6,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[6,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[6,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[6,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[6,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[6,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[6,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[6,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190926_N05m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N05m25.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[7,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[7,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists("THESIS_Results(New)/Result5/190926_P05m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P05m25.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[7,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[7,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[7,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[7,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[7,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[7,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[7,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[7,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190926_N05m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N05m50.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[8,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[8,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
if(file.exists("THESIS_Results(New)/Result5/190926_P05m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P05m50.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[8,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[8,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[8,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[8,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[8,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[8,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[8,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[8,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}



if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=40_nd=",21,".Rdata"))){
  load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=40_nd=",21,".Rdata"))}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[9,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[9,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
# if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#   load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
if(file.exists(paste0("THESIS_Results(New)/Result5/190823_P20m00.Rdata"))){
  load(file=paste0("THESIS_Results(New)/Result5/190823_P20m00.Rdata"))
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[9,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[9,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[9,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[9,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[9,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[9,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[9,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[9,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190802_N20m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190802_N20m10.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[10,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[10,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
# if(file.exists("THESIS_Results(New)/Result5/190802_P20m10.Rdata")){
#   load(file="THESIS_Results(New)/Result5/190802_P20m10.Rdata")
if(file.exists("THESIS_Results(New)/Result5/190823_P20m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190823_P20m10.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[10,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[10,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[10,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[10,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[10,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[10,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[10,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[10,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190802_N20m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190802_N20m25.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[11,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[11,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
# if(file.exists("THESIS_Results(New)/Result5/190802_P20m25.Rdata")){
#   load(file="THESIS_Results(New)/Result5/190802_P20m25.Rdata")
if(file.exists("THESIS_Results(New)/Result5/190823_P20m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190823_P20m25.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[11,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[11,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[11,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[11,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[11,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[11,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[11,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[11,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}
if(file.exists("THESIS_Results(New)/Result5/190802_N20m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190802_N20m50.Rdata")}
Pred.result = summary(mcmc.list(NSBSR.PredMat))
Y = simdata$Y
Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
plot(Y,Yhat)
Result5outMat[12,1] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
Result5outMat[12,6] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
# if(file.exists("THESIS_Results(New)/Result5/190802_P20m50.Rdata")){
#   load(file="THESIS_Results(New)/Result5/190802_P20m50.Rdata")
if(file.exists("THESIS_Results(New)/Result5/190823_P20m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190823_P20m50.Rdata")
  Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
  Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
  Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
  Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
  
  Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[12,2] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[12,7] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[12,3] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[12,8] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[12,4] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[12,9] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
  Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
  Result5outMat[12,5] <- sqrt(mean((Y-Yhat)^2,na.rm=T)) %>% round(3)
  Result5outMat[12,10] <- (lm(Y~Yhat) %>% summary)$adj.r.squared %>% round(2)
}

Result5outMat









