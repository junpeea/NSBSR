
rm(list = ls());
setwd("D:/NSBSR_FINAL")
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  


rm(list = ls());
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  

ND = 100; D = 16; grid.size = c(32,32)
result.cp.Mat = matrix(NA,nrow=ND,ncol=7*6)
result.cor.Mat = matrix(NA,nrow=ND,ncol=7*6)

# par(mfrow=c(2,2))
# fields::image.plot(matrix(Y,nrow=sqrt(length(Y))))
# fields::image.plot(matrix(Yhat,nrow=sqrt(length(Y))))
# fields::image.plot(matrix((Yhat.low < Y) * (Y < Yhat.upp),nrow=sqrt(length(Y))))
# plot(Y,Yhat)

for(nd in 1:100){
  ND = 100
  setwd("D:/NSBSR_FINAL")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result$quantiles[1:simdata$nloc,1]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result$quantiles[1:simdata$nloc,5]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    
    result.cp.Mat[nd,(0*6+1)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(0*6+1)] <- cor(Y,Yhat)
    
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result00$quantiles[1:simdata$nloc,1]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result00$quantiles[1:simdata$nloc,5]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(0*6+2)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(0*6+2)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result01$quantiles[1:simdata$nloc,1]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result01$quantiles[1:simdata$nloc,5]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(0*6+3)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(0*6+3)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result05$quantiles[1:simdata$nloc,1]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result05$quantiles[1:simdata$nloc,5]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(0*6+4)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(0*6+4)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result20$quantiles[1:simdata$nloc,1]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result20$quantiles[1:simdata$nloc,5]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(0*6+5)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(0*6+5)] <- cor(Y,Yhat)
    
  }
  
  cat(nd,"1/7 Nugget procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result$quantiles[1:simdata$nloc,1]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result$quantiles[1:simdata$nloc,5]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    
    result.cp.Mat[nd,(1*6+1)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(1*6+1)] <- cor(Y,Yhat)
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result00$quantiles[1:simdata$nloc,1]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result00$quantiles[1:simdata$nloc,5]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(1*6+2)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(1*6+2)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result01$quantiles[1:simdata$nloc,1]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result01$quantiles[1:simdata$nloc,5]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(1*6+3)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(1*6+3)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result05$quantiles[1:simdata$nloc,1]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result05$quantiles[1:simdata$nloc,5]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(1*6+4)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(1*6+4)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result20$quantiles[1:simdata$nloc,1]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result20$quantiles[1:simdata$nloc,5]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(1*6+5)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(1*6+5)] <- cor(Y,Yhat)
    
  }
  
  cat(nd,"2/7 Matern01 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result$quantiles[1:simdata$nloc,1]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result$quantiles[1:simdata$nloc,5]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    
    result.cp.Mat[nd,(2*6+1)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(2*6+1)] <- cor(Y,Yhat)
  }
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result00$quantiles[1:simdata$nloc,1]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result00$quantiles[1:simdata$nloc,5]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(2*6+2)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(2*6+2)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result01$quantiles[1:simdata$nloc,1]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result01$quantiles[1:simdata$nloc,5]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(2*6+3)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(2*6+3)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result05$quantiles[1:simdata$nloc,1]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result05$quantiles[1:simdata$nloc,5]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(2*6+4)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(2*6+4)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result20$quantiles[1:simdata$nloc,1]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result20$quantiles[1:simdata$nloc,5]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(2*6+5)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(2*6+5)] <- cor(Y,Yhat)
    
  }
  
  cat(nd,"3/7 Matern05 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result$quantiles[1:simdata$nloc,1]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result$quantiles[1:simdata$nloc,5]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    
    result.cp.Mat[nd,(3*6+1)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(3*6+1)] <- cor(Y,Yhat)
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result00$quantiles[1:simdata$nloc,1]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result00$quantiles[1:simdata$nloc,5]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(3*6+2)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(3*6+2)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result01$quantiles[1:simdata$nloc,1]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result01$quantiles[1:simdata$nloc,5]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(3*6+3)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(3*6+3)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result05$quantiles[1:simdata$nloc,1]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result05$quantiles[1:simdata$nloc,5]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(3*6+4)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(3*6+4)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result20$quantiles[1:simdata$nloc,1]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result20$quantiles[1:simdata$nloc,5]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(3*6+5)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(3*6+5)] <- cor(Y,Yhat)
  }
  
  cat(nd,"4/7 SqExp05 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result$quantiles[1:simdata$nloc,1]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result$quantiles[1:simdata$nloc,5]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    
    result.cp.Mat[nd,(4*6+1)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(4*6+1)] <- cor(Y,Yhat)
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result00$quantiles[1:simdata$nloc,1]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result00$quantiles[1:simdata$nloc,5]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(4*6+2)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(4*6+2)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result01$quantiles[1:simdata$nloc,1]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result01$quantiles[1:simdata$nloc,5]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(4*6+3)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(4*6+3)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result05$quantiles[1:simdata$nloc,1]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result05$quantiles[1:simdata$nloc,5]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(4*6+4)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(4*6+4)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result20$quantiles[1:simdata$nloc,1]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result20$quantiles[1:simdata$nloc,5]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(4*6+5)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(4*6+5)] <- cor(Y,Yhat)
  }
  
  cat(nd,"5/7 SqExp15 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result$quantiles[1:simdata$nloc,1]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result$quantiles[1:simdata$nloc,5]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    
    result.cp.Mat[nd,(5*6+1)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(5*6+1)] <- cor(Y,Yhat)
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result00$quantiles[1:simdata$nloc,1]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result00$quantiles[1:simdata$nloc,5]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(5*6+2)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(5*6+2)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result01$quantiles[1:simdata$nloc,1]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result01$quantiles[1:simdata$nloc,5]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(5*6+3)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(5*6+3)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result05$quantiles[1:simdata$nloc,1]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result05$quantiles[1:simdata$nloc,5]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(5*6+4)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(5*6+4)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result20$quantiles[1:simdata$nloc,1]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result20$quantiles[1:simdata$nloc,5]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(5*6+5)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(5*6+5)] <- cor(Y,Yhat)
  }
  
  cat(nd,"6/7 Matern20 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result$quantiles[1:simdata$nloc,1]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result$quantiles[1:simdata$nloc,5]+Pred.result$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    
    result.cp.Mat[nd,(6*6+1)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(6*6+1)] <- cor(Y,Yhat)
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
    Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
    Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
    Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result00$quantiles[1:simdata$nloc,1]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result00$quantiles[1:simdata$nloc,5]+Pred.result00$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(6*6+2)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(6*6+2)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result01$quantiles[1:simdata$nloc,1]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result01$quantiles[1:simdata$nloc,5]+Pred.result01$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(6*6+3)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(6*6+3)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result05$quantiles[1:simdata$nloc,1]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result05$quantiles[1:simdata$nloc,5]+Pred.result05$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(6*6+4)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(6*6+4)] <- cor(Y,Yhat)
    
    Y = simdata$Y
    Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.low = Pred.result20$quantiles[1:simdata$nloc,1]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),1]
    Yhat.upp = Pred.result20$quantiles[1:simdata$nloc,5]+Pred.result20$quantiles[(1+simdata$nloc):(2*simdata$nloc),5]
    result.cp.Mat[nd,(6*6+5)] <- mean((Yhat.low < Y) * (Y < Yhat.upp))
    result.cor.Mat[nd,(6*6+5)] <- cor(Y,Yhat)
  }
  
  cat(nd,"7/7 Gauss procedure complete!");print(Sys.time())
  
}
summary(result.cp.Mat)
summary(result.cor.Mat)


colMeans(result.cp.Mat,na.rm=T) %>% round(3)
colMeans(result.cor.Mat,na.rm=T) %>% round(3)

write.csv(result.cp.Mat,file=paste0(Sys.Date(),"result_cp_Mat.csv"))
write.csv(result.cor.Mat,file=paste0(Sys.Date(),"result_cor_Mat.csv"))



# rm(list = ls());
# source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
# source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  
# 
# ND = 100; D = 16; grid.size = c(32,32)
# mspe    = matrix(NA,nrow=ND,ncol=7*6)
# for(nd in 1:100){
#   
#   ND = 100
#   setwd("D:/NSBSR_FINAL")
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result = summary(mcmc.list(NSBSR.PredMat))
#     
#     Y = simdata$Y
#     Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#     mspe[nd,(0*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
#     Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
#     Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
#     Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
#     
#     Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(0*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(0*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(0*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(0*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   cat(nd,"1/7 Nugget procedure complete!",Sys.time(),"\n")
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result = summary(mcmc.list(NSBSR.PredMat))
#     
#     Y = simdata$Y
#     Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#     mspe[nd,(1*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
#     Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
#     Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
#     Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
#     
#     Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(1*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(1*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(1*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(1*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   cat(nd,"2/7 Matern01 procedure complete!",Sys.time(),"\n")
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result = summary(mcmc.list(NSBSR.PredMat))
#     
#     Y = simdata$Y
#     Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#     mspe[nd,(2*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))  
#   }
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
#     Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
#     Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
#     Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
#     
#     Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(2*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(2*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(2*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(2*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))  
#   }
#   
#   cat(nd,"3/7 Matern05 procedure complete!",Sys.time(),"\n")
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result = summary(mcmc.list(NSBSR.PredMat))
#     
#     Y = simdata$Y
#     Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#     mspe[nd,(3*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
#     Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
#     Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
#     Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
#     
#     Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(3*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(3*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(3*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(3*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   cat(nd,"4/7 SqExp05 procedure complete!",Sys.time(),"\n")
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result = summary(mcmc.list(NSBSR.PredMat))
#     
#     Y = simdata$Y
#     Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#     mspe[nd,(4*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
#     Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
#     Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
#     Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
#     
#     Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(4*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(4*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(4*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(4*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   cat(nd,"5/7 SqExp15 procedure complete!",Sys.time(),"\n")
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result = summary(mcmc.list(NSBSR.PredMat))
#     
#     Y = simdata$Y
#     Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#     mspe[nd,(5*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
#     Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
#     Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
#     Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
#     
#     Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(5*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(5*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(5*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(5*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   cat(nd,"6/7 Matern20 procedure complete!",Sys.time(),"\n")
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result = summary(mcmc.list(NSBSR.PredMat))
#     
#     Y = simdata$Y
#     Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#     mspe[nd,(6*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
#     load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
#     Pred.result00 = summary(mcmc.list(PBSR00.PredMat))
#     Pred.result01 = summary(mcmc.list(PBSR01.PredMat))
#     Pred.result05 = summary(mcmc.list(PBSR05.PredMat))
#     Pred.result20 = summary(mcmc.list(PBSR20.PredMat))
#     
#     Yhat = Pred.result00$statistics[1:simdata$nloc,1]+Pred.result00$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(6*6+2)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result01$statistics[1:simdata$nloc,1]+Pred.result01$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(6*6+3)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result05$statistics[1:simdata$nloc,1]+Pred.result05$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(6*6+4)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#     Yhat = Pred.result20$statistics[1:simdata$nloc,1]+Pred.result20$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
#     mspe[nd,(6*6+5)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
#   }
#   
#   cat(nd,"7/7 Gauss procedure complete!",Sys.time(),"\n")
#   
# }
# mspe
# boxplot(mspe)


