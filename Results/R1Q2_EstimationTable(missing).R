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
ComputeMSE = function(mcmc_list){
  Esti.beta = rbind(mcmc_list[[1]],mcmc_list[[2]],mcmc_list[[3]]); Esti.beta = Esti.beta[,1:3]
  True.beta = matrix(rep(c(0.01,0.02,0.03),nrow(Esti.beta)),ncol=3,byrow=TRUE)
  (Esti.beta-True.beta)^2 %>% sqrt %>% colMeans %>% as.numeric
}


ResultoutMat = matrix(NA,nrow=12,ncol=5*3)
nd <- 22; ND <- 100
truebeta = c(0.01,0.02,0.03)
if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
  load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[1:3,1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[1:3,1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
  
}

if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
  load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[1:3,2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[1:3,3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[1:3,4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[1:3,5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[1:3,2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[1:3,3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[1:3,4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[1:3,5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result5/190926_N01m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N01m10.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[4:6,1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[4:6,1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result5/190926_P01m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P01m10.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[4:6,2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[4:6,3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[4:6,4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[4:6,5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[4:6,2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[4:6,3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[4:6,4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[4:6,5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}


if(file.exists("THESIS_Results(New)/Result5/190926_N01m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N01m25.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[7:9,1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[7:9,1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result5/190926_P01m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P01m25.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[7:9,2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[7:9,3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[7:9,4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[7:9,5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[7:9,2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[7:9,3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[7:9,4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[7:9,5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result5/190926_N01m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N01m50.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[10:12,1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[10:12,1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result5/190926_P01m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P01m50.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[10:12,2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[10:12,3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[10:12,4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[10:12,5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[10:12,2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[10:12,3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[10:12,4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[10:12,5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=100_nd=21.Rdata")){
  load(file="THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=100_nd=21.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[1:3,5+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[1:3,5+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,5+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=100_nd=21.Rdata")){
  load(file="THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=100_nd=21.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[1:3,5+2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[1:3,5+3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[1:3,5+4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[1:3,5+5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[1:3,5+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,5+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,5+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,5+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,5+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[1:3,5+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[1:3,5+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[1:3,5+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190926_N05m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N05m10.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[4:6,5+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[4:6,5+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,5+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result5/190926_P05m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P05m10.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[4:6,5+2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[4:6,5+3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[4:6,5+4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[4:6,5+5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[4:6,5+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,5+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,5+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,5+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,5+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[4:6,5+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[4:6,5+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[4:6,5+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190926_N05m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N05m25.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[7:9,5+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[7:9,5+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,5+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190926_P05m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P05m25.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[7:9,5+2] <- Esti.result00$statistics[1:3,1] - truebeta
  ResultoutMat[7:9,5+3] <- Esti.result01$statistics[1:3,1] - truebeta
  ResultoutMat[7:9,5+4] <- Esti.result05$statistics[1:3,1] - truebeta
  ResultoutMat[7:9,5+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,5+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,5+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,5+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,5+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,5+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[7:9,5+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[7:9,5+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[7:9,5+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
  
}
if(file.exists("THESIS_Results(New)/Result5/190926_N05m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_N05m50.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[10:12,5+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[10:12,5+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,5+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190926_P05m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190926_P05m50.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[10:12,5+2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[10:12,5+3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[10:12,5+4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[10:12,5+5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[10:12,5+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,5+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,5+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,5+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,5+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[10:12,5+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[10:12,5+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[10:12,5+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}

if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=100_nd=",21,".Rdata"))){
  load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=100_nd=",21,".Rdata"))
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[1:3,10+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[1:3,10+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,10+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}
if(file.exists(paste0("THESIS_Results(New)/Result5/190823_P20m00.Rdata"))){
  load(file=paste0("THESIS_Results(New)/Result5/190823_P20m00.Rdata"))
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[1:3,10+2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[1:3,10+3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[1:3,10+4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[1:3,10+5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[1:3,10+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,10+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,10+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,10+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[1:3,10+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[1:3,10+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[1:3,10+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[1:3,10+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190802_N20m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190802_N20m10.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[4:6,10+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[4:6,10+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,10+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190823_P20m10.Rdata")){
  load(file="THESIS_Results(New)/Result5/190823_P20m10.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[4:6,10+2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[4:6,10+3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[4:6,10+4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[4:6,10+5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[4:6,10+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,10+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,10+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,10+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[4:6,10+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[4:6,10+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[4:6,10+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[4:6,10+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190802_N20m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190802_N20m25.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[7:9,10+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[7:9,10+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,10+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190823_P20m25.Rdata")){
  load(file="THESIS_Results(New)/Result5/190823_P20m25.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[7:9,10+2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[7:9,10+3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[7:9,10+4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[7:9,10+5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[7:9,10+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,10+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,10+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,10+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[7:9,10+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[7:9,10+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[7:9,10+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[7:9,10+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}
if(file.exists("THESIS_Results(New)/Result5/190802_N20m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190802_N20m50.Rdata")
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  ResultoutMat[10:12,10+1] <- Esti.result$statistics[1:3,1]
  # ResultoutMat[10:12,10+1] <- Esti.result$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,10+1] <- ComputeMSE(mcmc.list(NSBSR.EstiMat))
}

if(file.exists("THESIS_Results(New)/Result5/190823_P20m50.Rdata")){
  load(file="THESIS_Results(New)/Result5/190823_P20m50.Rdata")
  Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
  Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
  Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
  Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
  ResultoutMat[10:12,10+2] <- Esti.result00$statistics[1:3,1]
  ResultoutMat[10:12,10+3] <- Esti.result01$statistics[1:3,1]
  ResultoutMat[10:12,10+4] <- Esti.result05$statistics[1:3,1]
  ResultoutMat[10:12,10+5] <- Esti.result20$statistics[1:3,1]
  # ResultoutMat[10:12,10+2] <- Esti.result00$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,10+3] <- Esti.result01$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,10+4] <- Esti.result05$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,10+5] <- Esti.result20$statistics[1:3,1] - truebeta
  # ResultoutMat[10:12,10+2] <- ComputeMSE(mcmc.list(PBSR00.EstiMat))
  # ResultoutMat[10:12,10+3] <- ComputeMSE(mcmc.list(PBSR01.EstiMat))
  # ResultoutMat[10:12,10+4] <- ComputeMSE(mcmc.list(PBSR05.EstiMat))
  # ResultoutMat[10:12,10+5] <- ComputeMSE(mcmc.list(PBSR20.EstiMat))
}

ResultoutMat

write.csv(round(ResultoutMat,3),file=paste0(Sys.Date(),"ResultoutMat.csv"))






