########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : PSBSR

# Script name : PBSR_Main

# Date : 201811xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################

rm(list = ls())

####################################################################### IMPORT LIBRARIES & Functions ###
source("Library/NSBSR_ImportLibrary.R")
source("Library/NSBSR_Func.R")
source("Library/NSBSR_Func_MCMC_ver4.R")

################################################################################# INITIAL PROCEDURES ###
X.seed <- 1; Y.seed <- 1                               # Seed settings
grid.size = c(32,32); obs.delta = c(2,2); obs.prob=1   # Environment settings
ND = 30                                                # Experiment settings

source("Common/Get2dGridEnvironment.R")                # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                 # Generate 2-dimension lattice grid covariates
source("Common/SimulateGaussianProcess.R")             # Simulate 2-dimension Gaussian Process -> save dataset
# source("Suppl/SaveSampleImages.R")          # Describe & Save images of the simulated Y


##################################################################################### IMPORT DATASET ###
load(file=paste0("Dataset/Isotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata"))
# load(file=paste0("Dataset/AnIsotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))

############################################################################### ITERATION PROCEDURES ###

for(nd in 1:10){
  
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  
  ################################################################################## Initialize PBSR ###
  source("Common/Initialize_PBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=100,n.burn=0)
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
  
  start = Sys.time()
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
            PredResult <- Predict_PBSR(PBSR00.new.mcmcobj,simdata)
            PBSR00.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                            as.vector(PredResult$result2),
                                                            as.vector(PredResult$result3))

            PredResult <- Predict_PBSR(PBSR01.new.mcmcobj,simdata)
            # #######################################################################
            # # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
            # # png(paste0("temp_chain",chain,"iter",iter,".png"))
            # par(mfrow=c(2,2))
            # X = simdata$X; Y = simdata$Y
            # Yhat = PredResult$result1 + PredResult$result2
            # brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
            # score = sqrt(mean((Y-Yhat)^2))
            # 
            # Ymat =  matrix(Y,   ncol=sqrt(simdata$nloc))
            # Yhmat = matrix(Yhat,ncol=sqrt(length(Yhat)))
            # 
            # col <- rev(rainbow(77)[1:63])
            # 
            # fields::image.plot(Ymat,main="TRUE")
            # fields::image.plot(Yhmat,main=paste0("MSE=",round(score,2)))
            # # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
            # plot(Ymat,Yhmat);abline(b=1,lty=2,col=2)
            # par(mfrow=c(1,1))
            # fields::image.plot(PredResult$result4)
            
            # # dev.off()
            # #######################################################################
            PBSR01.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                            as.vector(PredResult$result2),
                                                            as.vector(PredResult$result3))

            PredResult <- Predict_PBSR(PBSR05.new.mcmcobj,simdata)
            # #######################################################################
            # # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
            # # png(paste0("temp_chain",chain,"iter",iter,".png"))
            # par(mfrow=c(2,2))
            # Y = data0$Y
            # X = cbind(data0$X0,data0$X1,data0$X2)
            # Yhat = Y
            # Yhat[which(data0$bdd2==0)] = PredResult$result1 + PredResult$result2
            # Yhat[which(data0$obs2==1)] = Y[which(data0$obs2==1)]
            # brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
            # score = sqrt(mean((Y-Yhat)^2))
            # 
            # Y = matrix(data0$Y, ncol=sqrt(length(data0$Y)))
            # Yhat = matrix(Yhat,ncol=sqrt(length(Yhat)))
            # 
            # col <- rev(rainbow(77)[1:63])
            # 
            # fields::image.plot(Y,main="TRUE")
            # fields::image.plot(Yhat,main="PRED(0.5)")
            # # fields::image.plot(Y,breaks=brk,col=col,main="PBSR")
            # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
            # plot(Y,Yhat);abline(b=1,lty=2,col=2)
            # # dev.off()
            # #######################################################################
            PBSR05.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                            as.vector(PredResult$result2),
                                                            as.vector(PredResult$result3))
            # 
            PredResult <- Predict_PBSR(PBSR20.new.mcmcobj,simdata)
            # #######################################################################
            # # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
            # # png(paste0("temp_chain",chain,"iter",iter,".png"))
            # par(mfrow=c(2,2))
            # Y = data0$Y
            # X = cbind(data0$X0,data0$X1,data0$X2)
            # Yhat = Y
            # Yhat[which(data0$bdd2==0)] = PredResult$result1 + PredResult$result2
            # Yhat[which(data0$obs2==1)] = Y[which(data0$obs2==1)]
            # brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
            # score = sqrt(mean((Y-Yhat)^2))
            # 
            # Y = matrix(data0$Y, ncol=sqrt(length(data0$Y)))
            # Yhat = matrix(Yhat,ncol=sqrt(length(Yhat)))
            # 
            # col <- rev(rainbow(77)[1:63])
            # 
            # fields::image.plot(Y,main="TRUE")
            # fields::image.plot(Yhat,main="PRED(2.0)")
            # # fields::image.plot(Y,breaks=brk,col=col,main="PBSR")
            # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
            # plot(Y,Yhat);abline(b=1,lty=2,col=2)
            # # dev.off()
            # #######################################################################
            PBSR20.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                            as.vector(PredResult$result2),
                                                            as.vector(PredResult$result3))


          }

          PBSR00.mcmcobj0 <- PBSR00.new.mcmcobj
          PBSR01.mcmcobj0 <- PBSR01.new.mcmcobj
          PBSR05.mcmcobj0 <- PBSR05.new.mcmcobj
          PBSR20.mcmcobj0 <- PBSR20.new.mcmcobj

          # cat("PBSR_chain",chain,"iter",iter," ",date(),'\n')

          cat("PBSR_chain",chain,"/",n.chain,"iter",iter,"/",n.iter,
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

  save.image(paste0("Result210907/",Sys.Date(),"_PBSR(nd=",nd,").RData"))

}

record = Sys.time() - start
record


