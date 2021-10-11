########################################################################################################
# Thesis : Bayesian Spatial Process Regression with Nonparametric modeling of spectral densities

# Project name : NSBSR_main

# Date : 201807xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################


############################################################################################## SETTINGS
rm(list = ls());
setwd("C:/Users/user/Desktop/180102_LCY_THESIS") # set your main directory
source(file="PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="NSBSR_Func.R")             # Functions for Nonarametric Spectral Bayesian Spatial Regression

######################################## Available versions
# source(file="NSBSR_Func_MCMC_ver1.R")      # Version1. Isotropic theta
source(file="NSBSR_Func_MCMC_ver2.R")        # Version2. common version
# source(file="NSBSR_Func_MCMC_ver2_Rcpp.R")   # Version2. common version
# source(file="NSBSR_Func_MCMC_ver3.R")      # Version3. 2-sided rho-hyper-par
########################################

D = 32;ND = 1
########################################################################################################

for(nd in 1:ND){
  
  # Generate 2-dimension lattice grid map for NSBSR
  map = Get2dgridmap(seed=1,lgs=D,resol=2)
  
  # Simulate Gaussian Random Field over the lattice map
  #####################################################################################################
  # aniso.option=RMangle(angle=0,ratio=1)############################################ ISOTROPY EXPERIMENT
  aniso.option=RMangle(angle=pi/4,ratio=0.5)######################################## ANISOTROPY EXPERIMENT
  # simdata = Simulate.Nugget(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  # simdata = Simulate.SqExp(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  # simdata = Simulate.SqExp(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(1.5,10),aniso.option=aniso.option)
  # simdata = Simulate.Matern(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.1,10),aniso.option=aniso.option)
  simdata = Simulate.Matern(map,seed=3,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  # simdata = Simulate.Matern(map,seed=1,nd=nd,sig2_eps=1,cor.par=c(2.0,10),aniso.option=aniso.option)
  # simdata = Simulate.Gauss(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  ###################################################################################################
  # Fundamental settings...
  data0     = simdata$DB;                            geodata0  = as.geodata(data0 ,coords.col=1:2,data.col=13,covar.col=9:11)
  data00    = simdata$DB[which(simdata$DB$obs==1),]; geodata00 = as.geodata(data00,coords.col=1:2,data.col=13,covar.col=9:11)
  # rho.vec  = exp(seq(log(0.01),log(10),length.out=100))
  # rho.vec  = exp(seq(log(0.1),log(1),length.out=100))
  # rho.vec  = exp(seq(log(0.0001),log(0.01),length.out=500))
  # rho.vec  = exp(seq(log(0.005),log(0.05),length.out=100))
  # rho.vec  = exp(seq(log(0.001),log(0.01),length.out=100))
  rho.vec  = exp(seq(log(0.0001),log(0.005),length.out=100))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  
  # TRUE
  loc = data00[,c("x","y")]-floor(D/2)
  # iso=TRUE
  iso=FALSE
  Rstar = matrix(c(-1/sqrt(5),1/sqrt(5),0,1),2,2)
  Rstar = matrix(c(-1/sqrt(5),0,1/sqrt(5),1),2,2)
  an.loc  = t(Rstar%*%t(loc))
  
  # truespd  <- matern.specdens(omega=loc,param=c((simdata$true.cov.pars[3]/D),simdata$true.cov.pars[1]))
  truespd  <- matern.specdens(omega=an.loc,param=c((simdata$true.cov.pars[3]/D),simdata$true.cov.pars[1]))
  iso.truespd <- truespd / sum(truespd)
  
  ################################################################################## Initialize NSBSR ###
  null.mcmcobj.list = SetMCMC(seed=200,nd=nd,n.chain=3,n.iter=10,n.burn=1)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = list(n.chain);NSBSR.PredMat = list(n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  ################################################################################## Initialize PBSR ###
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=10,n.burn=5)
  n.chain = length(null.mcmcobj.list)
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
  
  # ############################################### MCMC - PBSR #################################################
  for(chain in 1:n.chain){
# 
#     PBSR01.mcmcobj0 <- init.mcmcobj.PBSR01.list[[chain]]
#     PBSR05.mcmcobj0 <- init.mcmcobj.PBSR05.list[[chain]]
#     PBSR20.mcmcobj0 <- init.mcmcobj.PBSR20.list[[chain]]
#     n.iter = PBSR01.mcmcobj0$n.iter
#     n.burn = PBSR01.mcmcobj0$n.burn
# 
#     for(iter in 1:n.iter){
# 
#       PBSR01.new.mcmcobj <- Gibbs_PBSR(PBSR01.mcmcobj0,simdata)
#       PBSR05.new.mcmcobj <- Gibbs_PBSR(PBSR05.mcmcobj0,simdata)
#       PBSR20.new.mcmcobj <- Gibbs_PBSR(PBSR20.mcmcobj0,simdata)
# 
#       if(iter>n.burn){
# 
#         PBSR01.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR01.new.mcmcobj$Parset$beta,
#                                                         PBSR01.new.mcmcobj$Parset$kappa,
#                                                         PBSR01.new.mcmcobj$Parset$nu2,
#                                                         PBSR01.new.mcmcobj$Parset$sigma2,
#                                                         PBSR01.new.mcmcobj$Parset$range)
#         PBSR05.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR05.new.mcmcobj$Parset$beta,
#                                                         PBSR05.new.mcmcobj$Parset$kappa,
#                                                         PBSR05.new.mcmcobj$Parset$nu2,
#                                                         PBSR05.new.mcmcobj$Parset$sigma2,
#                                                         PBSR05.new.mcmcobj$Parset$range)
#         PBSR20.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR20.new.mcmcobj$Parset$beta,
#                                                         PBSR20.new.mcmcobj$Parset$kappa,
#                                                         PBSR20.new.mcmcobj$Parset$nu2,
#                                                         PBSR20.new.mcmcobj$Parset$sigma2,
#                                                         PBSR20.new.mcmcobj$Parset$range)
# 
#         PredResult <- Predict_PBSR(PBSR01.new.mcmcobj,simdata)
#         #######################################################################
#         # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
#         # png(paste0("temp_chain",chain,"iter",iter,".png"))
#         par(mfrow=c(2,2))
#         Y = data0$Y
#         X = cbind(data0$X0,data0$X1,data0$X2)
#         Yhat = Y
#         Yhat[which(data0$bdd2==0)] = PredResult$result1 + PredResult$result2
#         Yhat[which(data0$obs2==1)] = Y[which(data0$obs2==1)]
#         brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#         score = sqrt(mean((Y-Yhat)^2))
# 
#         Y = matrix(data0$Y, ncol=sqrt(length(data0$Y)))
#         Yhat = matrix(Yhat,ncol=sqrt(length(Yhat)))
# 
#         col <- rev(rainbow(77)[1:63])
# 
#         fields::image.plot(Y,main="TRUE")
#         fields::image.plot(Yhat,main="PRED(0.1)")
#         # fields::image.plot(Y,breaks=brk,col=col,main="PBSR")
#         fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
#         plot(Y,Yhat);abline(b=1,lty=2,col=2)
#         # dev.off()
#         #######################################################################
#         PBSR01.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
#                                                         as.vector(PredResult$result2),
#                                                         as.vector(PredResult$result3))
#         PredResult <- Predict_PBSR(PBSR05.new.mcmcobj,simdata)
#         #######################################################################
#         # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
#         # png(paste0("temp_chain",chain,"iter",iter,".png"))
#         par(mfrow=c(2,2))
#         Y = data0$Y
#         X = cbind(data0$X0,data0$X1,data0$X2)
#         Yhat = Y
#         Yhat[which(data0$bdd2==0)] = PredResult$result1 + PredResult$result2
#         Yhat[which(data0$obs2==1)] = Y[which(data0$obs2==1)]
#         brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#         score = sqrt(mean((Y-Yhat)^2))
# 
#         Y = matrix(data0$Y, ncol=sqrt(length(data0$Y)))
#         Yhat = matrix(Yhat,ncol=sqrt(length(Yhat)))
# 
#         col <- rev(rainbow(77)[1:63])
# 
#         fields::image.plot(Y,main="TRUE")
#         fields::image.plot(Yhat,main="PRED(0.5)")
#         # fields::image.plot(Y,breaks=brk,col=col,main="PBSR")
#         fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
#         plot(Y,Yhat);abline(b=1,lty=2,col=2)
#         # dev.off()
#         #######################################################################
#         PBSR05.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
#                                                         as.vector(PredResult$result2),
#                                                         as.vector(PredResult$result3))
#         PredResult <- Predict_PBSR(PBSR20.new.mcmcobj,simdata)
#         #######################################################################
#         # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
#         # png(paste0("temp_chain",chain,"iter",iter,".png"))
#         par(mfrow=c(2,2))
#         Y = data0$Y
#         X = cbind(data0$X0,data0$X1,data0$X2)
#         Yhat = Y
#         Yhat[which(data0$bdd2==0)] = PredResult$result1 + PredResult$result2
#         Yhat[which(data0$obs2==1)] = Y[which(data0$obs2==1)]
#         brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
#         score = sqrt(mean((Y-Yhat)^2))
# 
#         Y = matrix(data0$Y, ncol=sqrt(length(data0$Y)))
#         Yhat = matrix(Yhat,ncol=sqrt(length(Yhat)))
# 
#         col <- rev(rainbow(77)[1:63])
# 
#         fields::image.plot(Y,main="TRUE")
#         fields::image.plot(Yhat,main="PRED(2.0)")
#         # fields::image.plot(Y,breaks=brk,col=col,main="PBSR")
#         fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
#         plot(Y,Yhat);abline(b=1,lty=2,col=2)
#         # dev.off()
#         #######################################################################
#         PBSR20.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
#                                                         as.vector(PredResult$result2),
#                                                         as.vector(PredResult$result3))
# 
# 
# 
# 
#       }
# 
#       PBSR01.mcmcobj0 <- PBSR01.new.mcmcobj
#       PBSR05.mcmcobj0 <- PBSR05.new.mcmcobj
#       PBSR20.mcmcobj0 <- PBSR20.new.mcmcobj
# 
#       cat("PBSR_chain",chain,"iter",iter," ",date(),'\n')
# 
#     }
# 
#     PBSR01.EstiMat[[chain]] <- PBSR01.new.mcmcobj$EstiMat
#     PBSR01.PredMat[[chain]] <- PBSR01.new.mcmcobj$PredMat
#     PBSR05.EstiMat[[chain]] <- PBSR05.new.mcmcobj$EstiMat
#     PBSR05.PredMat[[chain]] <- PBSR05.new.mcmcobj$PredMat
#     PBSR20.EstiMat[[chain]] <- PBSR20.new.mcmcobj$EstiMat
#     PBSR20.PredMat[[chain]] <- PBSR20.new.mcmcobj$PredMat
#     
    ############################################ MCMC - NSBSR ###############################################

    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter
    n.burn = mcmcobj0$n.burn

    for(iter in 1:n.iter){

      # mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,simdata,a=1,b=1)
      mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,simdata,a=100,b=100)
      mcmcobj2 <- Update_Theta_NSBSR(mcmcobj1,simdata)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=0.01,d=0.01)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=0.01,d=1)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=1,d=1)
      mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=10,d=10)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=100,d=10)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=100,d=100)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=100,d=1000)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=1000,d=100)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=100,d=10000)
      mcmcobj4 <- Update_Mixture_NSBSR(mcmcobj3,simdata)

      new.mcmcobj <- mcmcobj4

      if(iter>n.burn){

        new.mcmcobj$EstiMat[(iter-n.burn),] <- c(new.mcmcobj$Parset$beta,
                                                 new.mcmcobj$Parset$tau.e,
                                                 new.mcmcobj$Parset$theta,
                                                 new.mcmcobj$Parset$rho.th)
        if((iter-n.burn)>1){
          Esti.summary = summary(mcmc(new.mcmcobj$EstiMat[1:(iter-n.burn),]))
          meantheta = Esti.summary$statistics[(simdata$n.beta+1+1):(simdata$n.beta+1+new.mcmcobj$n.theta),1]
          mat.theta  = matrix(c(meantheta,rev(meantheta)),ncol=sqrt(simdata$map$nobs))
          PredResult <- Predict_NSBSR(new.mcmcobj,mat.theta,simdata)

          #######################################################################
          # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
          # png(paste0("temp_chain",chain,"iter",iter,".png"))
          par(mfrow=c(2,3))
          Y = data0$Y
          X = cbind(data0$X0,data0$X1,data0$X2)
          Yhat = Y
          Yhat[which(data0$bdd2==0)] = PredResult$result1 + PredResult$result2
          Yhat[which(data0$obs2==1)] = Y[which(data0$obs2==1)]
          brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
          score = sqrt(mean((Y-Yhat)^2))

          Y = matrix(data0$Y, ncol=sqrt(length(data0$Y)))
          Yhat = matrix(Yhat,ncol=sqrt(length(Yhat)))

          col <- rev(rainbow(77)[1:63])

          fields::image.plot(Y,main="TRUE")
          fields::image.plot(Yhat,main=paste0("MSE=",round(score,2)))
          # fields::image.plot(Yhat,breaks=brk,col=col,main=paste0("MSE=",round(score,2)))

          plot(Y,Yhat);abline(b=1,lty=2,col=2)

          fields::image.plot(matrix(log(iso.truespd),ncol=sqrt(length(iso.truespd))),zlim=c(-14,0))

          if((iter-n.burn)>1){
            Esti.summary = summary(mcmc(new.mcmcobj$EstiMat[1:(iter-n.burn),]))
            meantheta = Esti.summary$statistics[(simdata$n.beta+1+1):(simdata$n.beta+1+new.mcmcobj$n.theta),1]
            comp1  = matrix(c(meantheta,rev(meantheta)),ncol=sqrt(simdata$map$nobs))
            rotate <- function(x) t(apply(x, 2, rev));
            comp2 <- comp1; if(iso==TRUE) comp2 <- rotate(comp1)
            # mat.theta = cbind(rbind(exp(comp1[5:8,5:8]),exp(comp2[1:4,5:8])),rbind(exp(comp2[5:8,1:4]),exp(comp1[1:4,1:4])))
            # mat.theta = cbind(rbind(exp(comp1[7:12,7:12]),exp(comp2[1:6,7:12])),rbind(exp(comp2[7:12,1:6]),exp(comp1[1:6,1:6])))
            # mat.theta = cbind(rbind(exp(comp1[9:16,9:16]),exp(comp1[1:8,9:16])),rbind(exp(comp1[9:16,1:8]),exp(comp1[1:8,1:8])))
            mat.theta = cbind(rbind(exp(comp1[(D/2+1):D,(D/2+1):D]),exp(comp2[1:(D/2),(D/2+1):D])),rbind(exp(comp2[(D/2+1):D,1:(D/2)]),exp(comp1[1:(D/2),1:(D/2)])))
            mat.theta = (mat.theta-max(mat.theta))*min(log(iso.truespd))/min((mat.theta-max(mat.theta)))
            fields::image.plot(mat.theta,main="SPD",zlim=c(-14,0))

            lambda = new.mcmcobj$Parset$lambda
            mat.lambda = matrix(lambda,ncol=D)/sum(lambda)
            # rotate <- function(x) t(apply(x, 2, rev))
            # mat.lambda <- rotate(mat.lambda)
            fields::image.plot(mat.lambda,main="SPD2")
            # dev.off()
            #######################################################################
          }

          new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                   as.vector(PredResult$result2),
                                                   as.vector(PredResult$result3))

        }






      }

      mcmcobj0 <- new.mcmcobj


      cat("NSBSR_chain",chain,"iter",iter," ",date(),'\n')

    }

    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat

  }
  
  # record = Sys.time() - start
  # setwd("Z:/개인 파일/NSBSR_JYB/Result/Result1")
  # # save.image(paste0("Iso_Data",nd,"_Mat(0.5)","_(iter=",n.iter,"_lgs=",map$gridsize,")",".Rdata"))
  # save.image(paste0("AnIso_Data",nd,"_Mat(0.5)","_(iter=",n.iter,"_lgs=",map$gridsize,")",".Rdata"))
  
}











