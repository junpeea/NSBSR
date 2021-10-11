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
setwd("C:/Users/user/Desktop/NSBSR_FINAL")
source("Library/NSBSR_ImportLibrary.R")
source("Library/NSBSR_Func.R")

################################################################################# INITIAL PROCEDURES ###
X.seed <- 1; Y.seed <- 1                               # Seed settings
grid.size = c(32,32); obs.delta = c(2,2)               # Environment settings
ND = 1                                      # Experiment settings
source("Common/Get2dGridEnvironment.R")               # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
source("Common/SimulateGaussianProcess.R")            # Simulate 2-dimension Gaussian Process -> save dataset
# source("Suppl/Save_Sample_Images.R")          # Describe & Save images of the simulated Y


############################################################################### ITERATION PROCEDURES ###

start = Sys.time()
for(nd in 1:ND){
  
  dir = paste0("Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd)
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  load(file=paste0("Dataset/Isotropic/",dir,".Rdata"))
  
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=E[[nd]], Y=Y[[nd]], name=filename)
  
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
  source("Library/NSBSR_Func_MCMC_ver5.R")
  
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
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result2")
  save.image(paste0("Iso",dir,".Rdata"))   
  
  
}

record = Sys.time() - start
record

for(nd in 1:ND){
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result2")
  load(file=paste0("IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))
  Pred.result = summary(mcmc.list(NSBSR.PredMat))
  
  Y = simdata$Y
  Yhat = Pred.result$statistics[1:simdata$nloc,1] + Pred.result$statistics[(simdata$nloc+1):(2*simdata$nloc),1]
  
  pdf(paste0(Sys.Date(),"True_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  fields::image.plot(matrix(Y,   sqrt(simdata$nloc),sqrt(simdata$nloc)),axes=F,cex.main=2.5,legend.mar=3.2,zlim=c(-2.1,2.1))
  dev.off()
  
  pdf(paste0(Sys.Date(),"Pred_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  fields::image.plot(matrix(Yhat,sqrt(simdata$nloc),sqrt(simdata$nloc)),axes=F,cex.main=2.5,legend.mar=3.2,zlim=c(-2.1,2.1))
  dev.off()
  
}




