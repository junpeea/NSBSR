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
setwd("D:/NSBSR_FINAL")
source("Library/NSBSR_ImportLibrary.R")
source("Library/NSBSR_Func.R")

################################################################################# INITIAL PROCEDURES ###
X.seed <- 1; Y.seed <- 1                               # Seed settings
grid.size = c(64,64); obs.delta = c(2,2)               # Environment settings
ND = 100 ; obs.prob=1.0                                     # Experiment settings
# ND = 30
mode="ISO"
source("Common/Get2dGridEnvironment.R")               # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
source("Common/SimulateGaussianProcess.R")            # Simulate 2-dimension Gaussian Process -> save dataset
# source("Suppl/Save_Sample_Images.R")          # Describe & Save images of the simulated Y


##################################################################################### IMPORT DATASET ###
# if(mode=="AnISO"){load(file=paste0("Dataset/AnIsotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata"))
#   }else{load(file=paste0("Dataset/Isotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata"))}
if(mode=="AnISO"){load(file=paste0("Dataset/AnIsotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata"))
}else{load(file=paste0("Dataset/Isotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata"))}

############################################################################### ITERATION PROCEDURES ###


start = Sys.time()
nd <- 1
for(nd in 1:1){
  
  # Load the Simulation dataset
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]),
                 X=X, E=El[[nd]], Y=Yl[[nd]], name=filename)
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=1000; d=100
  # rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.vec  = seq(0.1,10,length.out=10)
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  # source("Common/Initialize_NSBSR.R")
  # source("Common/Initialize_NSBSR_nutha.R")
  source("Common/Initialize_INBSR.R")
  # source("Common/Initialize_INBSR2.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=200,n.burn=9)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  # source("Library/NSBSR_Func_MCMC_ver8.R")
  source("Library/NSBSR_Func_MCMC_ver9.R")
  # source("Library/NSBSR_Func_MCMC_ver9(2).R")
  
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
      
      # meantheta = mcmcobj0$Parset$theta
      # n = sqrt(simdata$nobs)
      # if(n%%2==0){
      #   mat.theta  = matrix(c(meantheta,rev(meantheta)),ncol=n)
      # }else{
      #   mat.theta  = matrix(c(meantheta,rev(meantheta[1:(length(meantheta)-1)])),ncol=n)
      # }
      # nh = ceiling(n/2)
      # rearranged.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
      #                          cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
      # 
      # fields::image.plot(rearranged.theta,main=paste0("SPD(",iter,")"),axes=F)
      
      
      # cat("NSBSR_chain",chain,"iter",iter," ",date(),'\n')
      cat("NSBSR_chain",chain,"iter",iter,"tau.th",new.mcmcobj$Parset$tau.th,"rho.th",
          new.mcmcobj$Parset$rho.th1,new.mcmcobj$Parset$rho.th2,date(),'\n')
      
    }
    
    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat
    
  }
  
  if(mode=="AnISO"){save.image(file=paste0(Sys.Date(),"_AnIsotropic_Matern05_","Grid(",grid.size[1],",",grid.size[2],")_nd=",nd,".Rdata"))  
    }else{save.image(file=paste0(Sys.Date(),"_Isotropic_Matern05_","Grid(",grid.size[1],",",grid.size[2],")_nd=",nd,".Rdata"))  }

}

record = Sys.time() - start
record

library(coda)

nd <- 1
for(nd in 1:ND){
  # setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result1")
  setwd("D://NSBSR_FINAL/THESIS_Results(New)/Result1")
  load(file=paste0("190801_SPD_Iso32",".Rdata"))

  Esti.result = summary(mcmc.list(NSBSR.EstiMat))

  thetahat  = Esti.result$statistics[(3+1+1):(3+1+145)]
  mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
  nh = ceiling(n/2)
  re.mat.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
                           cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))

  pdf(paste0(Sys.Date(),"NSBSR-estimated_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  fields::image.plot(re.mat.theta,axes=F,zlim=c(-10,-1),cex.main=3.2,legend.mar=5.0)
  fields::image.plot(120*re.mat.theta+674,axes=F,cex.main=3.2,legend.mar=5.0)
  dev.off()


  n1=32; n2=32; alpha=0.5
  gp2=gp(c(n1,n2),matern.specdens,c(0.5,alpha))
  dens2=matern.specdens(gp2$omega,c(0.5,alpha),d=2)
  dens2mat = matrix(dens2,n1,n2)
  logdens2mat = log(dens2mat / sum(dens2mat))
  nh = ceiling(n/2)
  re.logdens2mat = rbind(cbind(logdens2mat[rev(1:nh),rev(1:nh)],logdens2mat[rev(1:nh),1:nh]),
                       cbind(logdens2mat[1:nh,rev(1:nh)],logdens2mat[1:nh,1:nh]))
  pdf(paste0(Sys.Date(),"TrueSPD_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  fields::image.plot(re.logdens2mat,axes=F,zlim=c(-10,-1),cex.main=3.2,legend.mar=5.0)
  dev.off()

}


nd <- 1
for(nd in 1:ND){
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result1")
  # load(file=paste0("190801_SPD_AnIso32",".Rdata"))
  load(file=paste0("190805_SPD_AnIso32",".Rdata"))

  Esti.result = summary(mcmc.list(NSBSR.EstiMat))

  thetahat  = Esti.result$statistics[(3+1+1):(3+1+145)]
  mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
  nh = ceiling(n/2)
  re.mat.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
                       cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
  pdf(paste0(Sys.Date(),"NSBSR-estimated_","AnIsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  # fields::image.plot(re.mat.theta,axes=F,cex.main=3.2,legend.mar=5.0)
  # fields::image.plot(160*re.mat.theta+899,axes=F,zlim=c(-12,-0.5),cex.main=3.2,legend.mar=5.0)
  fields::image.plot((re.mat.theta*180+1250),axes=F,cex.main=3.2,legend.mar=5.0)
  dev.off()

  library(spectralGP)
  library(dplyr)
  n1=32; n2=32; alpha=0.5
  gp2=gp(c(n1,n2),matern.specdens,c(0.5,alpha))
  A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)
  # A = diag(2)
  dens2=matern.specdens(as.matrix(gp2$omega)%*%A,c(0.5,alpha),d=2)
  dens2mat = matrix(dens2,n1,n2)
  logdens2mat = log(dens2mat / sum(dens2mat))
  nh = ceiling(n/2)
  re.logdens2mat = rbind(cbind(logdens2mat[rev(1:nh),rev(1:nh)],logdens2mat[rev(1:nh),1:nh]),
                         cbind(logdens2mat[1:nh,rev(1:nh)],logdens2mat[1:nh,1:nh]))
  pdf(paste0(Sys.Date(),"TrueSPD_","AnIsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  fields::image.plot(t(re.logdens2mat),axes=F,cex.main=3.2,legend.mar=5.0)
  # fields::image.plot(t(re.logdens2mat),zlim=c(-12,-0.5),axes=F,cex.main=3.2,legend.mar=5.0)
  dev.off()
}


nd <- 1
for(nd in 1:ND){
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result1")
  load(file=paste0("190801_SPD_Iso64",".Rdata"))

  Esti.result = summary(mcmc.list(NSBSR.EstiMat))

  thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
  mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
  nh = ceiling(n/2)
  re.mat.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
                       cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
  pdf(paste0(Sys.Date(),"NSBSR-estimated_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,"(nua).pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  fields::image.plot(re.mat.theta,axes=F,cex.main=3.2,legend.mar=5.0)
  # fields::image.plot(72*re.mat.theta+495,axes=F,zlim=c(-12,-0.5),cex.main=3.2,legend.mar=5.0)

  dev.off()

  library(spectralGP)
  library(dplyr)
  n1=64; n2=64; alpha=0.5
  gp2=gp(c(n1,n2),matern.specdens,c(0.25,alpha))
  A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)
  A = diag(2)
  dens2=matern.specdens(as.matrix(gp2$omega)%*%A,c(0.25,alpha),d=2)
  dens2mat = matrix(dens2,n1,n2)
  logdens2mat = log(dens2mat / sum(dens2mat))
  nh = ceiling(n/2)
  re.logdens2mat = rbind(cbind(logdens2mat[rev(1:nh),rev(1:nh)],logdens2mat[rev(1:nh),1:nh]),
                         cbind(logdens2mat[1:nh,rev(1:nh)],logdens2mat[1:nh,1:nh]))
  pdf(paste0(Sys.Date(),"TrueSPD_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  par(mfrow=c(1,1),mar=c(1,1,1,1))
  fields::image.plot(t(re.logdens2mat),axes=F,cex.main=3.2,legend.mar=5.0)
  # fields::image.plot(t(re.logdens2mat),zlim=c(-12,-0.5),axes=F,cex.main=3.2,legend.mar=5.0)
  dev.off()
}

# nd <- 1
# for(nd in 1:ND){
#   setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result1")
#   load(file=paste0("190801_SPD_AnIso64",".Rdata"))
# 
#   Esti.result = summary(mcmc.list(NSBSR.EstiMat))
# 
#   thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
#   mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
#   nh = ceiling(n/2)
#   re.mat.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
#                        cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
#   pdf(paste0(Sys.Date(),"NSBSR-estimated_","AnIsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,"(nua).pdf"),width=7*1,height=7*1)
#   par(mfrow=c(1,1),mar=c(1,1,1,1))
#   fields::image.plot(re.mat.theta,axes=F,cex.main=3.2,legend.mar=5.0)
#   # fields::image.plot(160*re.mat.theta+1109,axes=F,zlim=c(-12,-0.5),cex.main=3.2,legend.mar=5.0)
#   # fields::image.plot(63*re.mat.theta+431,axes=F,zlim=c(-12,-0.5),cex.main=3.2,legend.mar=5.0)
#   dev.off()
# 
#   library(spectralGP)
#   library(dplyr)
#   n1=32; n2=32; alpha=0.5
#   gp2=gp(c(n1,n2),matern.specdens,c(0.25,alpha))
#   A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)
#   # A = diag(2)
#   dens2=matern.specdens(as.matrix(gp2$omega)%*%A,c(0.25,alpha),d=2)
#   dens2mat = matrix(dens2,n1,n2)
#   logdens2mat = log(dens2mat / sum(dens2mat))
#   nh = ceiling(n/2)
#   re.logdens2mat = rbind(cbind(logdens2mat[rev(1:nh),rev(1:nh)],logdens2mat[rev(1:nh),1:nh]),
#                          cbind(logdens2mat[1:nh,rev(1:nh)],logdens2mat[1:nh,1:nh]))
#   pdf(paste0(Sys.Date(),"TrueSPD_","AnIsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
#   par(mfrow=c(1,1),mar=c(1,1,1,1))
#   fields::image.plot(t(re.logdens2mat),axes=F,cex.main=3.2,legend.mar=5.0)
#   # fields::image.plot(t(re.logdens2mat),zlim=c(-12,-0.5),axes=F,cex.main=3.2,legend.mar=5.0)
#   dev.off()
# }

#
# nd <- 1
# for(nd in 1:ND){
#   setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result1")
#   load(file=paste0("IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".Rdata"))
#
#   Esti.result = summary(mcmc.list(NSBSR.EstiMat))
#
#   thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
#   mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
#   nh = ceiling(n/2)
#   re.mat.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
#                        cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
#   pdf(paste0(Sys.Date(),"NSBSR-estimated_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
#   par(mfrow=c(1,1),mar=c(1,1,1,1))
#   # fields::image.plot(re.mat.theta,axes=F,cex.main=3.2,legend.mar=5.0)
#   fields::image.plot(160*re.mat.theta+1111,axes=F,zlim=c(-12,-0.5),cex.main=3.2,legend.mar=5.0)
#
#   dev.off()
#
#
#
#   library(spectralGP)
#   library(dplyr)
#   n1=32; n2=32; alpha=0.5
#   gp2=gp(c(n1,n2),matern.specdens,c(0.5,alpha))
#   # A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)
#   A = diag(2)
#   dens2=matern.specdens(as.matrix(gp2$omega)%*%A,c(0.5,alpha),d=2)
#   dens2mat = matrix(dens2,n1,n2)
#   logdens2mat = log(dens2mat / sum(dens2mat))
#   nh = ceiling(n/2)
#   re.logdens2mat = rbind(cbind(logdens2mat[rev(1:nh),rev(1:nh)],logdens2mat[rev(1:nh),1:nh]),
#                          cbind(logdens2mat[1:nh,rev(1:nh)],logdens2mat[1:nh,1:nh]))
#   pdf(paste0(Sys.Date(),"TrueSPD_","IsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
#   par(mfrow=c(1,1),mar=c(1,1,1,1))
#   # fields::image.plot(t(re.logdens2mat),axes=F,cex.main=3.2,legend.mar=5.0)
#   fields::image.plot(t(re.logdens2mat),zlim=c(-12,-0.5),axes=F,cex.main=3.2,legend.mar=5.0)
#   dev.off()
# }
