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
ND = 40; obs.prob=1
source("Common/Get2dGridEnvironment.R")               # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
source("Common/SimulateGaussianProcess.R")            # Simulate 2-dimension Gaussian Process -> save dataset
# source("Suppl/Save_Sample_Images.R")          # Describe & Save images of the simulated Y


dir = list(
  paste0("Gauss_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND),
  paste0("Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND),
  paste0("SqExp15_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND),
  paste0("SqExp05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND),
  paste0("Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND),
  paste0("Matern01_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND),
  paste0("Nugget_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND)
)

for(type in 1:7){
  # ##################################################################################### IMPORT DATASET (sample) ###
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  load(file=paste0("Dataset/Isotropic/",dir[[type]],".Rdata"))
  # load(file=paste0("Dataset/AnIsotropic/",dir[[type]],".Rdata"))
  
  ############################################################################### ITERATION PROCEDURES ###
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  # par(mfrow=c(1,1))
  # fields::image.plot(matrix(simdata$Y,30,30),col=colindex[6:68])
  
  ################################################################################## Initialize PBSR ###
  source("Common/Initialize_PBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=51,n.burn=1)
  n.chain = length(null.mcmcobj.list)
  PBSR00.EstiMat  = list(n.chain);PBSR00.PredMat  = list(n.chain)
  init.mcmcobj.PBSR00.list  = lapply( c(1:length(null.mcmcobj.list)),
                                      function(ind) Initialize_PBSR(null.mcmcobj.list[[ind]],simdata) )
  PBSR01.EstiMat  = list(n.chain);PBSR01.PredMat  = list(n.chain)
  init.mcmcobj.PBSR01.list  = lapply( c(1:length(null.mcmcobj.list)),
                                      function(ind) Initialize_PBSR(null.mcmcobj.list[[ind]],simdata,kappa=0.1) )
  PBSR05.EstiMat  = list(n.chain);PBSR05.PredMat  = list(n.chain)
  init.mcmcobj.PBSR05.list  = lapply( c(1:length(null.mcmcobj.list)),
                                      function(ind) Initialize_PBSR(null.mcmcobj.list[[ind]],simdata,kappa=0.5) )
  PBSR20.EstiMat  = list(n.chain);PBSR20.PredMat  = list(n.chain)
  init.mcmcobj.PBSR20.list  = lapply( c(1:length(null.mcmcobj.list)),
                                      function(ind) Initialize_PBSR(null.mcmcobj.list[[ind]],simdata,kappa=2.0) )
  
  ############################################# MCMC - PBSR #################################################
  source("Library/PBSR_Func_MCMC.R")
  ############################################# MCMC - PBSR #################################################
  chain <- 1
  for(chain in 1:n.chain){
    
    PBSR00.mcmcobj0 <- init.mcmcobj.PBSR00.list[[chain]]
    PBSR01.mcmcobj0 <- init.mcmcobj.PBSR01.list[[chain]]
    PBSR05.mcmcobj0 <- init.mcmcobj.PBSR05.list[[chain]]
    PBSR20.mcmcobj0 <- init.mcmcobj.PBSR20.list[[chain]]
    
    n.iter = PBSR00.mcmcobj0$n.iter
    n.burn = PBSR00.mcmcobj0$n.burn
    
    iter <- 1
    for(iter in 1:n.iter){
      
      PBSR00.new.mcmcobj <- Gibbs_PBSR(PBSR00.mcmcobj0,simdata,fix.kappa=FALSE)
      PBSR01.new.mcmcobj <- Gibbs_PBSR(PBSR01.mcmcobj0,simdata)
      PBSR05.new.mcmcobj <- Gibbs_PBSR(PBSR05.mcmcobj0,simdata)
      PBSR20.new.mcmcobj <- Gibbs_PBSR(PBSR20.mcmcobj0,simdata)
      
      if(iter>n.burn){
        
        PBSR00.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR00.new.mcmcobj$Parset$beta,
                                                        PBSR00.new.mcmcobj$Parset$kappa,
                                                        PBSR00.new.mcmcobj$Parset$nu2,
                                                        PBSR00.new.mcmcobj$Parset$sigma2,
                                                        PBSR00.new.mcmcobj$Parset$range)
        PBSR01.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR01.new.mcmcobj$Parset$beta,
                                                        PBSR01.new.mcmcobj$Parset$kappa,
                                                        PBSR01.new.mcmcobj$Parset$nu2,
                                                        PBSR01.new.mcmcobj$Parset$sigma2,
                                                        PBSR01.new.mcmcobj$Parset$range)
        PBSR05.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR05.new.mcmcobj$Parset$beta,
                                                        PBSR05.new.mcmcobj$Parset$kappa,
                                                        PBSR05.new.mcmcobj$Parset$nu2,
                                                        PBSR05.new.mcmcobj$Parset$sigma2,
                                                        PBSR05.new.mcmcobj$Parset$range)
        PBSR20.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR20.new.mcmcobj$Parset$beta,
                                                        PBSR20.new.mcmcobj$Parset$kappa,
                                                        PBSR20.new.mcmcobj$Parset$nu2,
                                                        PBSR20.new.mcmcobj$Parset$sigma2,
                                                        PBSR20.new.mcmcobj$Parset$range)
        
        # #######################################################################
        PredResult <- Predict_PBSR(PBSR00.new.mcmcobj,simdata)
        PBSR00.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                        as.vector(PredResult$result2),
                                                        as.vector(PredResult$result3))
        
        if(iter%%4==3){
          par(mfrow=c(2,2))
          X = simdata$X; Y = simdata$Y
          Yhat = PredResult$result1 + PredResult$result2
          brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
          score = sqrt(mean((Y-Yhat)^2))
          
          Ymat =  matrix(Y,   ncol=sqrt(simdata$nloc))
          Yhmat = matrix(Yhat,ncol=sqrt(length(Yhat)))
          
          col <- rev(rainbow(77)[1:63])
          
          fields::image.plot(Ymat,main="TRUE")
          fields::image.plot(Yhmat,main=paste0("PBSR00_","MSE=",round(score,2)))
          # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
          plot(Ymat,Yhmat);abline(b=1,lty=2,col=2)
          # par(mfrow=c(1,1))
          # fields::image.plot(PredResult$result4)
          
          # # dev.off()
          
        }
        
        
        
        # ####################################################################### 
        PredResult <- Predict_PBSR(PBSR01.new.mcmcobj,simdata)
        PBSR01.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                        as.vector(PredResult$result2),
                                                        as.vector(PredResult$result3))
        
        if(iter%%4==0){
          par(mfrow=c(2,2))
          X = simdata$X; Y = simdata$Y
          Yhat = PredResult$result1 + PredResult$result2
          brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
          score = sqrt(mean((Y-Yhat)^2))
          
          Ymat =  matrix(Y,   ncol=sqrt(simdata$nloc))
          Yhmat = matrix(Yhat,ncol=sqrt(length(Yhat)))
          
          col <- rev(rainbow(77)[1:63])
          
          fields::image.plot(Ymat,main="TRUE")
          fields::image.plot(Yhmat,main=paste0("PBSR01_","MSE=",round(score,2)))
          # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
          plot(Ymat,Yhmat);abline(b=1,lty=2,col=2)
          # par(mfrow=c(1,1))
          # fields::image.plot(PredResult$result4)
          
        }
        
        
        # #######################################################################
        PredResult <- Predict_PBSR(PBSR05.new.mcmcobj,simdata)
        PBSR05.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                        as.vector(PredResult$result2),
                                                        as.vector(PredResult$result3))
        
        if(iter%%4==1){
          par(mfrow=c(2,2))
          X = simdata$X; Y = simdata$Y
          Yhat = PredResult$result1 + PredResult$result2
          brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
          score = sqrt(mean((Y-Yhat)^2))
          
          Ymat =  matrix(Y,   ncol=sqrt(simdata$nloc))
          Yhmat = matrix(Yhat,ncol=sqrt(length(Yhat)))
          
          col <- rev(rainbow(77)[1:63])
          
          fields::image.plot(Ymat,main="TRUE")
          fields::image.plot(Yhmat,main=paste0("PBSR05_","MSE=",round(score,2)))
          # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
          plot(Ymat,Yhmat);abline(b=1,lty=2,col=2)
          # par(mfrow=c(1,1))
          # fields::image.plot(PredResult$result4)
        }
        
        
        # # dev.off()
        # #######################################################################
        PredResult <- Predict_PBSR(PBSR20.new.mcmcobj,simdata)
        PBSR20.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                        as.vector(PredResult$result2),
                                                        as.vector(PredResult$result3))
        
        if(iter%%4==2){
          par(mfrow=c(2,2))
          X = simdata$X; Y = simdata$Y
          Yhat = PredResult$result1 + PredResult$result2
          brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
          score = sqrt(mean((Y-Yhat)^2))
          
          Ymat =  matrix(Y,   ncol=sqrt(simdata$nloc))
          Yhmat = matrix(Yhat,ncol=sqrt(length(Yhat)))
          
          col <- rev(rainbow(77)[1:63])
          
          fields::image.plot(Ymat,main="TRUE")
          fields::image.plot(Yhmat,main=paste0("PBSR20_","MSE=",round(score,2)))
          # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
          plot(Ymat,Yhmat);abline(b=1,lty=2,col=2)
          # par(mfrow=c(1,1))
          # fields::image.plot(PredResult$result4)
        }
        
        
      }
      
      PBSR00.mcmcobj0 <- PBSR00.new.mcmcobj
      PBSR01.mcmcobj0 <- PBSR01.new.mcmcobj
      PBSR05.mcmcobj0 <- PBSR05.new.mcmcobj
      PBSR20.mcmcobj0 <- PBSR20.new.mcmcobj
      
      # cat("PBSR_chain",chain,"iter",iter," ",date(),'\n')
      
      cat(nd,"PBSR_chain",chain,"/",n.chain,"iter",iter,"/",n.iter,
          "alpha",PBSR00.new.mcmcobj$Parset$kappa," ",date(),'\n')
      
    }
    
    PBSR00.EstiMat[[chain]] <- PBSR00.new.mcmcobj$EstiMat
    PBSR00.PredMat[[chain]] <- PBSR00.new.mcmcobj$PredMat
    PBSR01.EstiMat[[chain]] <- PBSR01.new.mcmcobj$EstiMat
    PBSR01.PredMat[[chain]] <- PBSR01.new.mcmcobj$PredMat
    PBSR05.EstiMat[[chain]] <- PBSR05.new.mcmcobj$EstiMat
    PBSR05.PredMat[[chain]] <- PBSR05.new.mcmcobj$PredMat
    PBSR20.EstiMat[[chain]] <- PBSR20.new.mcmcobj$EstiMat
    PBSR20.PredMat[[chain]] <- PBSR20.new.mcmcobj$PredMat
  }
  
  # setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3")
  save.image(paste0("PBSR_",dir[[type]],"_nd=",nd,".Rdata"))  
  
}



