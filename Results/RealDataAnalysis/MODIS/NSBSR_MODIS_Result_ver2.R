########################################################################################################
########################################################################################################
# Thesis : Bayesian Spatial Process Regression with Nonparametric modeling of spectral densities
# Project name : NSBSR_MODIS_Result

# Date : 201808xx
# Author : YB JUN / Chae Young Lim
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
require(rgdal)
require(gdalUtils)
require(raster)
require(ggmap)


### LOAD MODIS dataset ##############################################################################################################
setwd("C:/Users/user/Desktop/180102_LCY_THESIS/MODIS_PROJECT")
sds <- list()
gdalinfo("MOD08_E3.A2017145.006.2017153180952.hdf"); sds[[1]] <- get_subdatasets("MOD08_E3.A2017145.006.2017153180952.hdf")
gdalinfo("MOD08_E3.A2017153.006.2017163121114.hdf"); sds[[2]] <- get_subdatasets("MOD08_E3.A2017153.006.2017163121114.hdf")
gdalinfo("MOD08_E3.A2017161.006.2017171120513.hdf"); sds[[3]] <- get_subdatasets("MOD08_E3.A2017161.006.2017171120513.hdf")
gdalinfo("MOD08_E3.A2017169.006.2017178120811.hdf"); sds[[4]] <- get_subdatasets("MOD08_E3.A2017169.006.2017178120811.hdf")
gdalinfo("MOD08_E3.A2017177.006.2017192120526.hdf"); sds[[5]] <- get_subdatasets("MOD08_E3.A2017177.006.2017192120526.hdf")
gdalinfo("MOD08_E3.A2017185.006.2017193082915.hdf"); sds[[6]] <- get_subdatasets("MOD08_E3.A2017185.006.2017193082915.hdf")
gdalinfo("MOD08_E3.A2017193.006.2017201083533.hdf"); sds[[7]] <- get_subdatasets("MOD08_E3.A2017193.006.2017201083533.hdf")
gdalinfo("MOD08_E3.A2017201.006.2017212170536.hdf"); sds[[8]] <- get_subdatasets("MOD08_E3.A2017201.006.2017212170536.hdf")
gdalinfo("MOD08_E3.A2017209.006.2017217085735.hdf"); sds[[9]] <- get_subdatasets("MOD08_E3.A2017209.006.2017217085735.hdf")
gdalinfo("MOD08_E3.A2017217.006.2017234125534.hdf"); sds[[10]] <- get_subdatasets("MOD08_E3.A2017217.006.2017234125534.hdf")
gdalinfo("MOD08_E3.A2017225.006.2017233084637.hdf"); sds[[11]] <- get_subdatasets("MOD08_E3.A2017225.006.2017233084637.hdf")
gdalinfo("MOD08_E3.A2017233.006.2017249151735.hdf"); sds[[12]] <- get_subdatasets("MOD08_E3.A2017233.006.2017249151735.hdf")
gdalinfo("MOD08_E3.A2017241.006.2017249084348.hdf"); sds[[13]] <- get_subdatasets("MOD08_E3.A2017241.006.2017249084348.hdf")
#####################################################################################################################################


mylist = list()
colindex = rainbow(77,alpha=0.5)
weekind <- 1
for(weekind in 1:13){
  
  # png(paste0("figure",weekind,".png"))
  # gdal_translate(sds[[weekind]][987], dst_dataset = paste0("Ozone_","Mean_WEEK",weekind,".tif"))
  # gdal_translate(sds[[dayind]][841], dst_dataset = paste0("Ozone_","sd_",(0903+dayind),".tif"))
  # gdal_translate(sds[[dayind]][842], dst_dataset = paste0("Ozone_","min_",(0903+dayind),".tif"))
  # gdal_translate(sds[[dayind]][843], dst_dataset = paste0("Ozone_","max_",(0903+dayind),".tif"))
  gdal_translate(sds[[weekind]][844], dst_dataset = paste0("Ozone_","QAmean_",weekind,".tif"))
  # gdal_translate(sds[[dayind]][845], dst_dataset = paste0("Ozone_","QAsd_",(0903+dayind),".tif"))
  # gdal_translate(sds[[dayind]][846], dst_dataset = paste0("Ozone_","HistCounts_",(0903+dayind),".tif"))
  # gdal_translate(sds[[dayind]][847], dst_dataset = paste0("Ozone_","ConfidenceHist_",(0903+dayind),".tif"))
  # rast <- raster(paste0("Ozone_","Mean_WEEK",weekind,".tif"))
  rast <- raster(paste0("Ozone_","QAmean_",weekind,".tif"))
  totdata <- as.matrix(rast);dim(totdata)
  mylist[[weekind]] <- log(totdata[(24+90):(53+90),(112+180):(141+180)]+1e-10)
  # sum(is.na(log(totdata)))
  # mylist[[weekind]] <- log(totdata+1e-10)
  # dim(mylist[[weekind]])
  # fields::image.plot(t(mylist[[weekind]]),main=paste0(2017,"week",weekind),cex.main=1.5,col=colindex)
  # dev.off()
  
}

# dim(rast)
# # fields::image.plot(rast)
# 
# plot.raster(rast,xlim=c((24+90),(53+90),ylim=c((112+180),(141+180))))


# image plot over googlemap

mymap  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "satellite");plot(mymap)
mymap2  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "terrain-labels");plot(mymap2)
mymap3 = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "toner-lines");plot(mymap3)


x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="lon",ylab="lat",main="")
fields::image.plot(x=x,y=seq(24,53,by=1),z=t(mylist[[1]]),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

#####################################################################################################################################
setwd("C:/Users/user/Desktop/180102_LCY_THESIS")
source(file="PBSR_Func_MCMC.R")
# source(file="NSBSR_Func_MCMC_ver1.R")
source(file="NSBSR_Func_MCMC_ver2.R")
source(file="NSBSR_Func.R")


Realdata = mylist[[3]]; Realdata[is.na(Realdata)] <- mean(Realdata,na.rm=TRUE)
par(mfrow=c(1,1))
x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="lon",ylab="lat",main="")
image(x=x,y=seq(24,53,by=1),z=t(Realdata),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

# Figure06
map = Get2dgridmap(seed=1,lgs=nrow(mylist[[1]])/2,resol=2)
  map$fieldgrid$X1 <- map$fieldgrid$x; map$fieldgrid$X2 <- map$fieldgrid$y
simdata = Simulate.Matern(map,seed=1,sig2_eps=1,nd=1,cor.par=c(0.5,10))
  simdata$DB$Y[which(simdata$DB$bdd2==0)] <- as.vector(Realdata[2:30,2:30])
  simdata$DB$Y[which(simdata$DB$bdd2==1)] <- NA

# par(mfrow=c(1,1))
# x = seq(111,141,by=1);y = 28/28*(x-113)+25
# plot(x,y,type="n",xlab="lon",ylab="lat",main="")
# image(x=x,y=seq(23,53,by=1),z=t(matrix(simdata$DB$Y,ncol=31)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
# par(new=T);plot(mymap2)
# par(new=T);plot(mymap3)


#Figure07  
map = Get2dgridmap(seed=1,lgs=nrow(mylist[[1]]),resol=2)
map$fieldgrid$X1 <- map$fieldgrid$x; map$fieldgrid$X2 <- map$fieldgrid$y
simdata = Simulate.Matern(map,seed=1,sig2_eps=1,nd=1,cor.par=c(0.5,10))
Y = rep(NA,nrow(simdata$DB)); Y[which(simdata$DB$obs2==1)]<-as.vector(Realdata); simdata$DB$Y <- Y

# par(mfrow=c(1,2))
# x = seq(112,141,by=1);y = 28/28*(x-113)+25
# plot(x,y,type="n",xlab="lon",ylab="lat",main="")
# fields::image.plot(x=x,y=seq(24,53,by=1),z=t(Realdata),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
# par(new=T);plot(mymap2)
# par(new=T);plot(mymap3)
# 
# temp = matrix(simdata$DB$Y,ncol=sqrt(length(Y)))
# temp = temp[3:nrow(temp),3:ncol(temp)]
# Y = temp
# 
# x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
# plot(x,y,type="n",xlab="lon",ylab="lat",main="")
# fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(Y),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
# par(new=T);plot(mymap2)
# par(new=T);plot(mymap3)


ND = 1
########################################################################################################

for(nd in 1:ND){
  
  # # Generate 2-dimension lattice grid map for NSBSR
  # map = Get2dgridmap(seed=1,lgs=D,resol=2)
  # 
  # # Simulate Gaussian Random Field over the lattice map
  # #####################################################################################################
  # aniso.option=RMangle(angle=0,ratio=1)############################################ ISOTROPY EXPERIMENT
  # # aniso.option=RMangle(angle=pi/4,ratio=2)######################################## ANISOTROPY EXPERIMENT
  # # simdata = Simulate.Nugget(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  # # simdata = Simulate.SqExp(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  # # simdata = Simulate.SqExp(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(1.5,10),aniso.option=aniso.option)
  # # simdata = Simulate.Matern(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.1,10),aniso.option=aniso.option)
  # simdata = Simulate.Matern(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  # # simdata = Simulate.Matern(map,seed=1,nd=nd,sig2_eps=1,cor.par=c(2.0,10),aniso.option=aniso.option)
  # # simdata = Simulate.Gauss(map,seed=1,sig2_eps=1,nd=nd,cor.par=c(0.5,10),aniso.option=aniso.option)
  # ###################################################################################################
  # Fundamental settings...
  data0     = simdata$DB;                            geodata0  = as.geodata(data0 ,coords.col=1:2,data.col=13,covar.col=9:11)
  data00    = simdata$DB[which(simdata$DB$obs==1),]; geodata00 = as.geodata(data00,coords.col=1:2,data.col=13,covar.col=9:11)
  # rho.vec  = exp(seq(log(0.0001),log(10),length.out=1000))
  # rho.vec  = exp(seq(log(0.1),log(1),length.out=100))
  # rho.vec  = exp(seq(log(0.01),log(0.5),length.out=500))
  # rho.vec  = exp(seq(log(0.005),log(0.05),length.out=100))
  # rho.vec  = exp(seq(log(0.0001),log(0.001),length.out=100))
  rho.vec  = exp(seq(log(0.0001),log(0.005),length.out=500))
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
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=10000,n.burn=9000)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = list(n.chain);NSBSR.PredMat = list(n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  ################################################################################## Initialize PBSR ###
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=10,n.burn=5)
  n.chain = length(null.mcmcobj.list)
  PBSR.EstiMat  = list(n.chain);PBSR.PredMat  = list(n.chain)
  init.mcmcobj.PBSR.list  = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_PBSR(null.mcmcobj.list[[ind]],simdata) )
  start = Sys.time()
  
  ############################################### MCMC - PBSR #################################################
  for(chain in 1:n.chain){
    
    PBSR.mcmcobj0 <- init.mcmcobj.PBSR.list[[chain]]
    n.iter = PBSR.mcmcobj0$n.iter
    n.burn = PBSR.mcmcobj0$n.burn
    
    for(iter in 1:n.iter){
      
      PBSR.new.mcmcobj <- Gibbs_PBSR(PBSR.mcmcobj0,simdata)
      
      if(iter>n.burn){
        
        PBSR.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR.new.mcmcobj$Parset$beta,
                                                      PBSR.new.mcmcobj$Parset$kappa,
                                                      PBSR.new.mcmcobj$Parset$nu2,
                                                      PBSR.new.mcmcobj$Parset$sigma2,
                                                      PBSR.new.mcmcobj$Parset$range)
        
        PredResult <- Predict_PBSR(PBSR.new.mcmcobj,simdata)
        
        #######################################################################
        # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
        # png(paste0("temp_chain",chain,"iter",iter,".png"))
        par(mfrow=c(2,2))
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
        fields::image.plot(Yhat,main="PRED")
        # fields::image.plot(Y,breaks=brk,col=col,main="PBSR")
        fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
        plot(Y,Yhat);abline(b=1,lty=2,col=2)
        # dev.off()
        #######################################################################
        
        PBSR.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                      as.vector(PredResult$result2),
                                                      as.vector(PredResult$result3))
        
      }
      
      PBSR.mcmcobj0 <- PBSR.new.mcmcobj
      
      cat("PBSR_chain",chain,"iter",iter," ",date(),'\n')
      
    }
    
    PBSR.EstiMat[[chain]] <- PBSR.new.mcmcobj$EstiMat
    PBSR.PredMat[[chain]] <- PBSR.new.mcmcobj$PredMat
    
    ############################################ MCMC - NSBSR ###############################################
    
    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter
    n.burn = mcmcobj0$n.burn
    
    for(iter in 1:n.iter){
      
      # mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,simdata,a=1,b=1)
      mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,simdata,a=100,b=25)
      mcmcobj2 <- Update_Theta_NSBSR(mcmcobj1,simdata)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=0.01,d=0.01)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=0.01,d=1)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=1,d=1)
      # mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=10,d=10)
      mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=100,d=100)
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
          mat.theta  = matrix(c(meantheta,rev(meantheta[2:length(meantheta)])),ncol=sqrt(simdata$map$nobs))
          # mat.theta  = matrix(c(meantheta,rev(meantheta)),ncol=sqrt(simdata$map$nobs))
          PredResult <- Predict_NSBSR(new.mcmcobj,mat.theta,simdata)
        }
        
        
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
          comp1  = matrix(c(meantheta,rev(meantheta[2:length(meantheta)])),ncol=sqrt(simdata$map$nobs))
          # comp1  = matrix(c(meantheta,rev(meantheta)),ncol=sqrt(simdata$map$nobs))
          rotate <- function(x) t(apply(x, 2, rev));
          comp2 <- comp1; if(iso==TRUE) comp2 <- rotate(comp1)
          mat.theta = cbind(rbind(exp(comp1[9:15,9:15]),exp(comp2[1:7,9:15])),rbind(exp(comp2[9:15,1:7]),exp(comp1[1:7,1:7])))
          # mat.theta = cbind(rbind(exp(comp1[(D/2+1):D,(D/2+1):D]),exp(comp2[1:(D/2),(D/2+1):D])),rbind(exp(comp2[(D/2+1):D,1:(D/2)]),exp(comp1[1:(D/2),1:(D/2)])))
          mat.theta = (mat.theta-max(mat.theta))*min(log(iso.truespd))/min((mat.theta-max(mat.theta)))
          fields::image.plot(mat.theta,main="SPD",zlim=c(-14,0))
          
          lambda = new.mcmcobj$Parset$lambda
          mat.lambda = matrix(lambda,ncol=15)/sum(lambda)
          # rotate <- function(x) t(apply(x, 2, rev))
          # mat.lambda <- rotate(mat.lambda)
          fields::image.plot(mat.lambda,main="SPD2")
        } 
        
        
        
        # dev.off()
        #######################################################################
        
        new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                 as.vector(PredResult$result2),
                                                 as.vector(PredResult$result3))
        
      }
      
      mcmcobj0 <- new.mcmcobj
      
      
      cat("NSBSR_chain",chain,"iter",iter," ",date(),'\n')
      
    }
    
    NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
    NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat
    
  }
  
  # Esti.summary = summary(mcmc.list(NSBSR.EstiMat))
  # Pred.summary = summary(mcmc.list(NSBSR.PredMat))
  # pdf("NSBSR_temp.pdf")
  # plot(mcmc.list(NSBSR.EstiMat))
  # dev.off()
  # PBSR.Esti.summary = summary(mcmc.list(PBSR.EstiMat))
  # PBSR.Pred.summary = summary(mcmc.list(PBSR.PredMat))
  # pdf("PBSR_temp.pdf")
  # plot(mcmc.list(PBSR.EstiMat))
  # dev.off()
  
  record = Sys.time() - start
  setwd("C:/Users/user/Desktop/180102_LCY_THESIS/180719_Experiment")
  save.image(paste0("Iso_Data",nd,"_Mat(0.5)","_(iter=",n.iter,"_lgs=",map$gridsize,")",".Rdata"))
  # save.image(paste0("AnIso_Data",nd,"_Mat(0.5)","_(iter=",n.iter,"_lgs=",map$gridsize,")",".Rdata"))
  
}






png(paste0(Sys.Date(),"MODIS_FIG1.png"),width=360*2,height=360*2)

par(mfrow=c(2,2))
colindex = rainbow(77,alpha=0.5)
brk = (seq(-4.5,4,length.out=78))

x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[1]])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

x = seq(112,141,by=2);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="TRAINING DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=2),z=(t(mylist[[1]][seq(1,30,by=2),seq(1,30,by=2)])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


Pred.summary = summary(mcmc.list(NSBSR.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2


# x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
x = seq(112,141,by=1);y = 28/28*(x-113)+25
# x = seq(113,141,by=1);y = 28/28*(x-113)+25
sc = mean(t(mylist[[1]])[2:30,2:30] - t(temp3),na.rm=T)^2
r2 = summary(lm(as.vector(t(mylist[[1]])[2:30,2:30]) ~ as.vector(t(temp3))))$r.squared
# plot(x,y,type="n",xlab="",ylab="",main=paste0("ESTIMATED (NSBSR) OZONE (",round(r2,2), ")(",round(sc,4),")"),xaxt='n',yaxt='n')
plot(x,y,type="n",xlab="",ylab="",main=paste0("ESTIMATED (NSBSR) OZONE"),xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
# fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(24,53,by=0.5),z=exp(t(temp3[2:31,2:31])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
fields::image.plot(x=x,y=seq(25,53,by=1),z=(t(temp3)),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


# Pred.summary = summary(mcmc.list(PBSR.PredMat))

# n = nrow(Pred.summary$statistics)/3
# temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
# temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
# temp3 = temp + temp2

# x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
x = seq(112,141,by=1);y = 28/28*(x-113)+25
# x = seq(113,141,by=1);y = 28/28*(x-113)+25
sc = mean(t(mylist[[1]])[2:30,2:30] - t(temp3),na.rm=T)^2
r2 = summary(lm(as.vector(t(mylist[[1]])[2:30,2:30]) ~ as.vector(t(temp3))))$r.squared
# plot(x,y,type="n",xlab="",ylab="",main=paste0("ESTIMATED (PBSR) OZONE (",round(r2,2), ")(",round(sc,4),")"),xaxt='n',yaxt='n')
plot(x,y,type="n",xlab="",ylab="",main=paste0("ESTIMATED (PBSR) OZONE"),xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
# fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(24,53,by=0.5),z=exp(t(temp3[2:31,2:31])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
fields::image.plot(x=x,y=seq(25,53,by=1),z=(t(temp3)),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

dev.off()



png(paste0(Sys.Date(),"MODIS_FIG2.png"),width=360*2,height=360*2)

par(mfrow=c(2,2))
colindex = rainbow(77,alpha=0.5)
brk = (seq(-4.5,4,length.out=78))

x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[1]])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


# Pred.summary = summary(mcmc.list(NSBSR.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,3],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),3],ncol=sqrt(n))
temp3 = temp + temp2



x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
# x = seq(112,141,by=1);y = 28/28*(x-113)+25
# x = seq(113,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")

# fields::image.plot(x=x,y=seq(24,53,by=0.5),z=exp(t(temp3[2:31,2:31])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(25,53,by=1),z=(t(temp3)),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2


x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
# x = seq(112,141,by=1);y = 28/28*(x-113)+25
# x = seq(113,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE (2.5% quantile)",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(24,53,by=1),z=exp(t(temp3[2:31,2:31])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(25,53,by=1),z=exp(t(temp3)),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,5],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),5],ncol=sqrt(n))
temp3 = temp + temp2
# fields::image.plot(temp)
# fields::image.plot(t(temp2))

x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
# x = seq(112,141,by=1);y = 28/28*(x-113)+25
# x = seq(113,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE (97.5% quantile)",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(24,53,by=1),z=exp(t(temp3[2:31,2:31])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(25,53,by=1),z=exp(t(temp3)),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)



dev.off()

