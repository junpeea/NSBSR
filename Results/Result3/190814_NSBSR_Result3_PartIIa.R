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
ND = 100; obs.prob=1.0
source("Common/Get2dGridEnvironment.R")               # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
source("Common/SimulateGaussianProcess.R")            # Simulate 2-dimension Gaussian Process -> save dataset
# source("Suppl/Save_Sample_Images.R")          # Describe & Save images of the simulated Y


# ##################################################################################### IMPORT DATASET (sample) ###
dir = paste0("Gauss_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND)
setwd("C:/Users/user/Desktop/NSBSR_FINAL")
load(file=paste0("Dataset/AnIsotropic/",dir,".Rdata"))

############################################################################### ITERATION PROCEDURES ###

for(nd in 1:10){
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=1000; d=10
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  source("Common/Initialize_NSBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=500,n.burn=450)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  source("Library/NSBSR_Func_MCMC_ver7.R")
  
  chain <- 1
  for(chain in 1:n.chain){
    
    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter;  n.burn = mcmcobj0$n.burn
    
    
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
      cat(nd,"NSBSR_chain",chain,"/",n.chain,"iter",iter,"/",n.iter,"tau.th",new.mcmcobj$Parset$tau.th,"rho.th",
          new.mcmcobj$Parset$rho.th1,new.mcmcobj$Parset$rho.th2,date(),'\n')
      
    }
    
    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat
    
  }
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3")
  save.image(paste0("NSBSR_",dir,"_nd=",nd,".Rdata"))   
  
}



# ##################################################################################### IMPORT DATASET (sample) ###
dir = paste0("Matern01_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND)
setwd("C:/Users/user/Desktop/NSBSR_FINAL")
load(file=paste0("Dataset/AnIsotropic/",dir,".Rdata"))

############################################################################### ITERATION PROCEDURES ###

for(nd in 1:10){
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=1000; d=10
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  source("Common/Initialize_NSBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=500,n.burn=450)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  source("Library/NSBSR_Func_MCMC_ver7.R")
  
  chain <- 1
  for(chain in 1:n.chain){
    
    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter;  n.burn = mcmcobj0$n.burn
    
    
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
      cat(nd,"NSBSR_chain",chain,"/",n.chain,"iter",iter,"/",n.iter,"tau.th",new.mcmcobj$Parset$tau.th,"rho.th",
          new.mcmcobj$Parset$rho.th1,new.mcmcobj$Parset$rho.th2,date(),'\n')
      
    }
    
    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat
    
  }
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3")
  save.image(paste0("NSBSR_",dir,"_nd=",nd,".Rdata"))   
  
}



# ##################################################################################### IMPORT DATASET (sample) ###
dir = paste0("Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND)
setwd("C:/Users/user/Desktop/NSBSR_FINAL")
load(file=paste0("Dataset/AnIsotropic/",dir,".Rdata"))

############################################################################### ITERATION PROCEDURES ###

for(nd in 1:10){
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=1000; d=10
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  source("Common/Initialize_NSBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=500,n.burn=450)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  source("Library/NSBSR_Func_MCMC_ver7.R")
  
  chain <- 1
  for(chain in 1:n.chain){
    
    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter;  n.burn = mcmcobj0$n.burn
    
    
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
      cat(nd,"NSBSR_chain",chain,"/",n.chain,"iter",iter,"/",n.iter,"tau.th",new.mcmcobj$Parset$tau.th,"rho.th",
          new.mcmcobj$Parset$rho.th1,new.mcmcobj$Parset$rho.th2,date(),'\n')
      
    }
    
    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat
    
  }
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3")
  save.image(paste0("NSBSR_",dir,"_nd=",nd,".Rdata"))   
  
}


# ##################################################################################### IMPORT DATASET (sample) ###
dir = paste0("Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND)
setwd("C:/Users/user/Desktop/NSBSR_FINAL")
load(file=paste0("Dataset/AnIsotropic/",dir,".Rdata"))

############################################################################### ITERATION PROCEDURES ###

for(nd in 1:10){
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=1000; d=10
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  source("Common/Initialize_NSBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=500,n.burn=450)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  source("Library/NSBSR_Func_MCMC_ver7.R")
  
  chain <- 1
  for(chain in 1:n.chain){
    
    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter;  n.burn = mcmcobj0$n.burn
    
    
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
      cat(nd,"NSBSR_chain",chain,"/",n.chain,"iter",iter,"/",n.iter,"tau.th",new.mcmcobj$Parset$tau.th,"rho.th",
          new.mcmcobj$Parset$rho.th1,new.mcmcobj$Parset$rho.th2,date(),'\n')
      
    }
    
    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat
    
  }
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result3")
  save.image(paste0("NSBSR_",dir,"_nd=",nd,".Rdata"))   
  
}

