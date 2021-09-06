########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Script name : NSBSR_Result3

# Date : 201906xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################


rm(list = ls());
setwd("C:/Users/user/Desktop/NSBSR_FINAL")
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  

ND = 100; D = 16; grid.size = c(32,32)
mspe    = matrix(NA,nrow=ND,ncol=7*6)

for(nd in 1:30){
  
  setwd("D:/noname")

  if(file.exists(paste0("Isotropic/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(0*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("Isotropic/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))

    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  cat(nd,"1/7 Nugget procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Isotropic/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(1*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("Isotropic/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))

    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  cat(nd,"2/7 Matern01 procedure complete!");print(Sys.time())

  if(file.exists(paste0("Isotropic/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(2*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  if(file.exists(paste0("Isotropic/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  cat(nd,"3/7 Matern05 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Isotropic/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(3*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("Isotropic/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))

    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  cat(nd,"4/7 SqExp05 procedure complete!");print(Sys.time())

  if(file.exists(paste0("Isotropic/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(4*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("Isotropic/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))

    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  cat(nd,"5/7 SqExp15 procedure complete!");print(Sys.time())

  if(file.exists(paste0("Isotropic/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(5*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("Isotropic/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))

    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  cat(nd,"6/7 Matern20 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Isotropic/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(6*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("Isotropic/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Isotropic/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))

    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  cat(nd,"7/7 Gauss procedure complete!");print(Sys.time())
  
}
mspe
boxplot(mspe)
summary(mspe)

ND = 100; D = 16; grid.size = c(32,32)
mspe    = matrix(NA,nrow=ND,ncol=7*6)
for(nd in 1:30){
  
  setwd("D:/noname")
  
  if(file.exists(paste0("Anisotropic/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(0*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("Anisotropic/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(0*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  cat(nd,"1/7 Nugget procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Anisotropic/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(1*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("Anisotropic/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(1*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  cat(nd,"2/7 Matern01 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Anisotropic/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(2*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))  
  }
  
  if(file.exists(paste0("Anisotropic/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(2*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))  
  }
  
  cat(nd,"3/7 Matern05 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Anisotropic/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(3*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("Anisotropic/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(3*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  cat(nd,"4/7 SqExp05 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Anisotropic/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(4*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("Anisotropic/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(4*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  cat(nd,"5/7 SqExp15 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Anisotropic/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(5*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("Anisotropic/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(5*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  cat(nd,"6/7 Matern20 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("Anisotropic/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(6*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("Anisotropic/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("Anisotropic/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    mspe[nd,(6*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  cat(nd,"7/7 Gauss procedure complete!");print(Sys.time())
  
}
mspe
boxplot(mspe)
summary(mspe)

# temp2 <- mspe
# temp2[1:40,seq(1,42,by=6)] <- mspe[,seq(1,42,by=6)]
# boxplot(temp2)
# mspe <- temp2

setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3")
pdf(paste0("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3/",Sys.Date(),"Figure3_1(Iso).pdf"),width=7*2,height=7*1)
par(mfrow=c(1,1),mar=c(3,3,1,1))
boxplot(sqrt(mspe[,-c(42:48)]),xaxt="n",col=rep(c(2,"lightblue",4,"darkblue","green",NA),7),
        ylim=c(0,1.2),xlab="True corr",ylab="RMSPE",main="",cex.main=3.2)
axis(side=1,at=seq(3,39,by=6),c("Indep","Matern(0,1)","Exponential","Pow.Exp(0.5)","Pow.Exp(1.5)","Matern(2,0)","Gauss"),font=2.0,cex.axis=1.5)
abline(h=0.4,lty=2);abline(h=0.8,lty=2);abline(h=1.2,lty=2);abline(h=1.6,lty=2);abline(h=2.0,lty=2)
for(i in 1:6){abline(v=(6*i),lty=1)}
legend("bottomleft",col=c(2,"lightblue",4,"darkblue","green"),c("NSBSR","PBSR(a=0.1)","PBSR(a=0.5)","PBSR(a=2.0)","PBSR(a unfixed)"),pch=rep(15,5),bg="white",cex=1.5)
dev.off()

setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3")
pdf(paste0("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3/",Sys.Date(),"Figure3_2(AnIso).pdf"),width=7*2,height=7*1)
par(mfrow=c(1,1),mar=c(3,3,1,1))
boxplot(sqrt(mspe[,-c(42:48)]),xaxt="n",col=rep(c(2,"lightblue",4,"darkblue","green",NA),7),
        ylim=c(0,1.2),xlab="True corr",ylab="RMSPE",main="",cex.main=3.2)
axis(side=1,at=seq(3,39,by=6),c("Indep","Matern(0,1)","Exponential","Pow.Exp(0.5)","Pow.Exp(1.5)","Matern(2,0)","Gauss"),font=2.0,cex.axis=1.5)
abline(h=0.4,lty=2);abline(h=0.8,lty=2);abline(h=1.2,lty=2);abline(h=1.6,lty=2);abline(h=2.0,lty=2)
for(i in 1:6){abline(v=(6*i),lty=1)}
legend("bottomleft",col=c(2,"lightblue",4,"darkblue","green"),c("NSBSR","PBSR(a=0.1)","PBSR(a=0.5)","PBSR(a=2.0)","PBSR(a unfixed)"),pch=rep(15,5),bg="white",cex=1.5)
dev.off()

