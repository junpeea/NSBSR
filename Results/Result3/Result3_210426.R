########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Script name : NSBSR_Result3

# Date : 201906xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################

rm(list = ls())
####################################################################### IMPORT LIBRARIES & Functions ###
source("Library/NSBSR_ImportLibrary.R")
source("Library/NSBSR_Func.R")

################################################################################# INITIAL PROCEDURES ###
X.seed <- 1; Y.seed <- 1                               # Seed settings
grid.size = c(32,32); obs.delta = c(2,2)               # Environment settings
source("Common/Get2dGridEnvironment.R")               # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
for(ND in 1:30){
  source("Common/SimulateGaussianProcess.R")            # Simulate 2-dimension Gaussian Process -> save dataset
}
# source("Suppl/Save_Sample_Images.R")          # Describe & Save images of the simulated Y


##################################################################################### IMPORT DATASET (sample) ###
nd <- 1
load(file=paste0("Dataset/Isotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))

############################################################################### ITERATION PROCEDURES ###

start = Sys.time()
for(nd in 1:ND){
  
  dir = paste0("Gauss_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd)
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  load(file=paste0("Dataset/Isotropic/",dir,".Rdata"))
  
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  # Set NSBSR Hyperparameters
  a=100; b=1; c=10000; d=1
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  source("Common/Initialize_NSBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=1000,n.burn=900)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  # source("Library/NSBSR_Func_MCMC_ver5.R")
  source("Library/NSBSR_Func_MCMC_ver7.R")
  
  chain <- 1
  for(chain in 1:n.chain){
    
    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter;  n.burn = mcmcobj0$n.burn
    
    iter <- 1
    for(iter in 1:n.iter){
      
      if(iter%%25==1 | iter>n.burn){
        mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,simdata,a=a,b=b)
      }else{mcmcobj1 <- mcmcobj0}
      mcmcobj2 <- Update_Theta_NSBSR(mcmcobj1,simdata)
      if(iter%%5==1){
        mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=c,d=d)
      }else{mcmcobj3 <- mcmcobj2}
      mcmcobj4 <- Update_Mixture_NSBSR(mcmcobj3,simdata)
      
      new.mcmcobj <- mcmcobj4
      
      if(iter>n.burn){
        
        new.mcmcobj$EstiMat[(iter-n.burn),] <- c(new.mcmcobj$Parset$beta,
                                                 new.mcmcobj$Parset$tau.e,
                                                 new.mcmcobj$Parset$theta,
                                                 new.mcmcobj$Parset$rho.th1,
                                                 new.mcmcobj$Parset$rho.th2)
        if((iter-n.burn)>1){
          
          Esti.summary = summary(mcmc(new.mcmcobj$EstiMat[1:(iter-n.burn),]))
          
          PredResult <- Predict_NSBSR(new.mcmcobj,Esti.summary,simdata)
          
          source("Suppl/Show_MCMC_NSBSR_Images.R")
          
          new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                   as.vector(PredResult$result2),
                                                   as.vector(PredResult$result3))
        }
        
        
      }
      
      mcmcobj0 <- new.mcmcobj
      
      
      # cat("NSBSR_chain",chain,"iter",iter," ",date(),'\n')
      cat("NSBSR_chain",chain,"iter",iter,"tau.th",new.mcmcobj$Parset$tau.th,"rho.th",
          new.mcmcobj$Parset$rho.th1,new.mcmcobj$Parset$rho.th2,date(),'\n')
      
    }
    
    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat
    
  }
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)")
  save.image(paste0(dir,".Rdata"))   
  
  
}

record = Sys.time() - start
record


rm(list = ls());
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  

ND = 100; D = 16; grid.size = c(32,32)
mspe    = matrix(NA,nrow=ND,ncol=7*6)

for(nd in 1:100){
  ND = 100
  setwd("D:/NSBSR_FINAL")

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(0*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(1*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(2*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(3*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(4*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(5*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))

    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(6*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }

  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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


rm(list = ls());
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  

ND = 100; D = 16; grid.size = c(32,32)
mspe    = matrix(NA,nrow=ND,ncol=7*6)
for(nd in 1:100){
  
  ND = 100
  setwd("D:/NSBSR_FINAL")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(0*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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
  
  cat(nd,"1/7 Nugget procedure complete!",Sys.time(),"\n")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(1*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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
  
  cat(nd,"2/7 Matern01 procedure complete!",Sys.time(),"\n")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(2*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))  
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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
  
  cat(nd,"3/7 Matern05 procedure complete!",Sys.time(),"\n")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(3*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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
  
  cat(nd,"4/7 SqExp05 procedure complete!",Sys.time(),"\n")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(4*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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
  
  cat(nd,"5/7 SqExp15 procedure complete!",Sys.time(),"\n")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(5*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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
  
  cat(nd,"6/7 Matern20 procedure complete!",Sys.time(),"\n")
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Pred.result = summary(mcmc.list(NSBSR.PredMat))
    
    Y = simdata$Y
    Yhat = Pred.result$statistics[1:simdata$nloc,1]+Pred.result$statistics[(1+simdata$nloc):(2*simdata$nloc),1]
    brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
    mspe[nd,(6*6+1)] <- sqrt(mean((Y-Yhat)^2,na.rm=T))
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ANISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
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
  
  cat(nd,"7/7 Gauss procedure complete!",Sys.time(),"\n")
  
}
mspe
boxplot(mspe)


pdf(paste0("D:/NSBSR_FINAL/THESIS_Results(New)/Result3/",Sys.Date(),"Figure3_1(Iso).pdf"),width=7*2,height=7*1)
par(mfrow=c(1,1),mar=c(3,3,1,1))
boxplot(sqrt(mspe[,-c(42:48)]),xaxt="n",col=rep(c(2,"lightblue",4,"darkblue","green",NA),7),
        ylim=c(0,1.2),xlab="True corr",ylab="RMSPE",main="",cex.main=3.2)
axis(side=1,at=seq(3,39,by=6),c("Indep","Matern(0,1)","Exponential","Pow.Exp(0.5)","Pow.Exp(1.5)","Matern(2,0)","Gauss"),font=2.0,cex.axis=1.5)
abline(h=0.4,lty=2);abline(h=0.8,lty=2);abline(h=1.2,lty=2);abline(h=1.6,lty=2);abline(h=2.0,lty=2)
for(i in 1:6){abline(v=(6*i),lty=1)}
legend("bottomleft",col=c(2,"lightblue",4,"darkblue","green"),c("NSBSR","PBSR(a unfixed)","PBSR(a=0.1)","PBSR(a=0.5)","PBSR(a=2.0)"),pch=rep(15,5),bg="white",cex=1.5)
dev.off()



library(ggplot2)
library(ggthemes)
temp = sqrt(mspe[,-c(42:48)])
temp2 = data.frame(cbind(as.vector(temp),sort(rep(c(1:41),100))))
temp3 = data.frame(temp2); colnames(temp3) <- c("rmspe","gp")

pdf(paste0("D:/NSBSR_FINAL/THESIS_Results(New)/Result3/",Sys.Date(),"Figure3_1(Iso).pdf"),width=7*2,height=7*1)

colind = rainbow(5)
temp5 = as.vector(matrix(rep(c(colind[c(4,3,5,2,1)],NA),100*7),nrow=100,byrow=TRUE))[1:4100]
ggplot(data=temp3, aes(x=gp,y=rmspe,group=gp))+ theme_hc() +
  theme(text=element_text(size=20)) +
  scale_x_continuous("",breaks=seq(3,39,by=6),
                     labels=c("Indep","Matern(0,1)","Exponential","Pow.Exp(0.5)","Pow.Exp(1.5)","Matern(2,0)","Gauss")) +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  geom_vline(xintercept=6, linetype="dashed", color = "black") +
  geom_vline(xintercept=12, linetype="dashed", color = "black") +
  geom_vline(xintercept=18, linetype="dashed", color = "black") +
  geom_vline(xintercept=24, linetype="dashed", color = "black") +
  geom_vline(xintercept=30, linetype="dashed", color = "black") +
  geom_vline(xintercept=36, linetype="dashed", color = "black") +
  geom_vline(xintercept=42, linetype="dashed", color = "black") +
  # geom_boxplot(aes(fill=c(c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"))))
  geom_boxplot(aes(fill=temp5)) +
  scale_fill_manual("",values=colind,
                    label=c("NSBSR","PBSR(a unfixed)","PBSR(a=0.1)","PBSR(a=0.5)","PBSR(a=2.0)")
  )

dev.off()



pdf(paste0("D:/NSBSR_FINAL/THESIS_Results(New)/Result3/",Sys.Date(),"Figure3_2(AnIso).pdf"),width=7*2,height=7*1)
par(mfrow=c(1,1),mar=c(3,3,1,1))
boxplot(sqrt(mspe[,-c(42:48)]),xaxt="n",col=rep(c(2,"lightblue",4,"darkblue","green",NA),7),
        ylim=c(0,1.5),xlab="True corr",ylab="RMSPE",main="",cex.main=3.2)
axis(side=1,at=seq(3,39,by=6),c("Indep","Matern(0,1)","Exponential","Pow.Exp(0.5)","Pow.Exp(1.5)","Matern(2,0)","Gauss"),font=2.0,cex.axis=1.5)
abline(h=0.4,lty=2);abline(h=0.8,lty=2);abline(h=1.2,lty=2);abline(h=1.6,lty=2);abline(h=2.0,lty=2)
for(i in 1:6){abline(v=(6*i),lty=1)}
legend("bottomleft",col=c(2,"lightblue",4,"darkblue","green"),c("NSBSR","PBSR(a unfixed)","PBSR(a=0.1)","PBSR(a=0.5)","PBSR(a=2.0)"),pch=rep(15,5),bg="white",cex=1.5)
dev.off()




pdf(paste0("D:/NSBSR_FINAL/THESIS_Results(New)/Result3/",Sys.Date(),"Figure3_2(AnIso).pdf"),width=7*2,height=7*1)

library(ggplot2)
library(ggthemes)
temp = sqrt(mspe[,-c(42:48)])
temp2 = data.frame(cbind(as.vector(temp),sort(rep(c(1:41),100))))
temp3 = data.frame(temp2); colnames(temp3) <- c("mspe","gp")

colind = rainbow(5)
temp5 = as.vector(matrix(rep(c(colind[c(4,3,5,2,1)],NA),100*7),nrow=100,byrow=TRUE))[1:4100]
ggplot(data=temp3, aes(x=gp,y=mspe,group=gp)) + theme_hc() +
  theme(text=element_text(size=20)) +
  scale_x_continuous("",breaks=seq(3,39,by=6),
                     labels=c("Indep","Matern(0,1)","Exponential","Pow.Exp(0.5)","Pow.Exp(1.5)","Matern(2,0)","Gauss")) +
  geom_vline(xintercept=0, linetype="dashed", color = "black") +
  geom_vline(xintercept=6, linetype="dashed", color = "black") +
  geom_vline(xintercept=12, linetype="dashed", color = "black") +
  geom_vline(xintercept=18, linetype="dashed", color = "black") +
  geom_vline(xintercept=24, linetype="dashed", color = "black") +
  geom_vline(xintercept=30, linetype="dashed", color = "black") +
  geom_vline(xintercept=36, linetype="dashed", color = "black") +
  geom_vline(xintercept=42, linetype="dashed", color = "black") +
  # geom_boxplot(aes(fill=c(c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"),
  #                     c("red","lightblue","blue","darkblue","green"))))
  geom_boxplot(aes(fill=temp5)) +
  scale_fill_manual("",values=colind,
                    label=c("NSBSR","PBSR(a unfixed)","PBSR(a=0.1)","PBSR(a=0.5)","PBSR(a=2.0)")
                    )

dev.off()

