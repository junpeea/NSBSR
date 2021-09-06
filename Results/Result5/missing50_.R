########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Script name : NSBSR_Main

# Date : 201811xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################

rm(list = ls())
####################################################################### IMPORT LIBRARIES & Functions ###
source("Library/NSBSR_ImportLibrary.R")
source("Library/NSBSR_Func.R")

################################################################################# INITIAL PROCEDURES ###
X.seed <- 1; Y.seed <- 1                              # Seed settings
grid.size = c(32,32); obs.delta = c(2,2)              # Environment settings
obs.prob = 0.5; ND = 33                               # Experiment settings

par(mfrow=c(2,3))
source("Common/Get2dGridEnvironment.R")               # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
source("Common/SimulateGaussianProcess.R")            # Simulate 2-dimension Gaussian Process -> save dataset
# source("Suppl/Save_Sample_Images.R")          # Describe & Save images of the simulated Y

#################################################################################### IMPORT DATASET ###
# load(file=paste0("Dataset/Isotropic/Matern01_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))
load(file=paste0("Dataset/Isotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata"))
# load(file=paste0("Dataset/AnIsotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))
# load(file=paste0("Dataset/Isotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))
# load(file=paste0("Dataset/AnIsotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))

############################################################################### ITERATION PROCEDURES ###

start = Sys.time()
nd <- 1
for(nd in 1:1){
  
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  # table(simdata$map$obs)
  
  # Set NSBSR Hyperparameters
  a=100; b=1; c=10000; d=1
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  # source("Common/Initialize_NSBSR.R")
  source("Common/Initialize_INBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=100,n.burn=10)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  source("Library/NSBSR_Func_MCMC_ver9.R")
  
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
      mcmcobj5 <- Update_Imputed.Y_NSBSR(mcmcobj4,simdata)
      
      # new.mcmcobj <- mcmcobj4
      new.mcmcobj <- mcmcobj5
      
      # par(mfrow=c(1,2))
      # # fields::image.plot(matrix(simdata$Y[which(simdata$map$obs>0)],ncol=17),zlim=c(-2,3.5))
      # # fields::image.plot(matrix(new.mcmcobj$sampleY,ncol=17),zlim=c(-2,3.5))
      # plot(new.mcmcobj$Parset$beta,type='o',ylim=c(-0.1,0.5))
      # fields::image.plot(matrix(c(new.mcmcobj$Parset$theta,rev(new.mcmcobj$Parset$theta[2:145])),ncol=17))
      
      
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
  
  
  
}

record = Sys.time() - start
record


