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
require(rgdal)
require(gdalUtils)
require(raster)
require(ggmap)


### LOAD MODIS dataset ##############################################################################################################
setwd("C:/Users/yjun/Desktop/WORK2021/180102_LCY_THESIS/MODIS_PROJECT")
load("C:/Users/yjun/Desktop/WORK2021/180102_LCY_THESIS/MODIS_PROJECT/180217_tempsave.RData")

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


# image plot over googlemap
mymap  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "satellite");plot(mymap)
mymap2  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "terrain-labels");plot(mymap2)
mymap3 = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "toner-lines");plot(mymap3)


Realdata = mylist[[3]]; Realdata[is.na(Realdata)] <- mean(Realdata,na.rm=TRUE)
par(mfrow=c(1,1))
x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="lon",ylab="lat",main="")
image(x=x,y=seq(24,53,by=1),z=t(Realdata),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

#####################################################################################################################################
# Estimation Procedure
setwd("C:/Users/yjun/Desktop/WORK2021/180102_LCY_THESIS")
source(file="PBSR_Func_MCMC.R")
# source(file="NSBSR_Func_MCMC_ver1.R")
source(file="NSBSR_Func_MCMC_ver2.R")
source(file="NSBSR_Func.R")

map = Get2dgridmap(seed=1,lgs=nrow(mylist[[1]])/2,resol=2)
map$fieldgrid$X1 <- map$fieldgrid$x; map$fieldgrid$X2 <- map$fieldgrid$y
simdata = Simulate.Matern(map,seed=1,sig2_eps=1,nd=1,cor.par=c(0.5,10))
simdata$DB$Y[which(simdata$DB$bdd2==0)] <- as.vector(Realdata[2:30,2:30])
simdata$DB$Y[which(simdata$DB$bdd2==1)] <- NA

par(mfrow=c(1,1))
x = seq(111,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="lon",ylab="lat",main="")
image(x=x,y=seq(23,53,by=1),z=t(matrix(simdata$DB$Y,ncol=31)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


########################################################################################################
setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702")
source("Library/NSBSR_ImportLibrary.R")
source("Library/NSBSR_Func.R")

################################################################################# INITIAL PROCEDURES ###
X.seed <- 1; Y.seed <- 1                                   # Seed settings
grid.size = c(32,32); obs.delta = c(2,2); obs.prob=1   # Environment settings
ND = 1                                                   # Experiment settings
source("Common/Get2dGridEnvironment.R")               # Generate 2-dimension lattice grid map
source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
source("Common/SimulateGaussianProcess.R")            # Simulate 2-dimension Gaussian Process -> save dataset

nd = 1
for(nd in 1:ND){
  
  # Load the Simulation dataset
  Yl = matrix(NA,31,31); Yl[1:30,1:30] <- Realdata
  mysimdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nrow(simdata$DB), nobs=sum(!is.na(simdata$DB$obs)), map=simdata$DB[,c(1,2,3,4,7,5)],
                 beta=beta, true.cov.par=c(cor.par[1],sig2_eps,cor.par[2]), X=simdata$DB[,c("X0","x","y")], E=rep(NA,nrow(simdata$DB)), Y=as.numeric(Yl), name=filename)
  colnames(mysimdata$X) <- c("X0","X1","X2")
  mysimdata$map$obs[which(mysimdata$map$obs==1)] <- 4
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=100; d=100
  rho.vec  = exp(seq(log(1e-03),log(1e+02),length.out=25))
  rho.grid = Supp_grid_search_rho(mysimdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  source("Common/Initialize_NSBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=100,n.burn=0)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],mysimdata) )
  
  ############################################ MCMC - NSBSR #################################################
  source("Library/NSBSR_Func_MCMC_ver14.R")
  
  chain <- 1
  for(chain in 1:n.chain){
    
    mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
    n.iter = mcmcobj0$n.iter;  n.burn = mcmcobj0$n.burn
    
    iter <- 1
    for(iter in 1:n.iter){
      
      # if(iter%%25==1 | iter>n.burn){
      if(1){
        mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,mysimdata,a=a,b=b)
      }else{mcmcobj1 <- mcmcobj0}
      mcmcobj2 <- Update_Theta_NSBSR(mcmcobj1,mysimdata)
      if(1){
        mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,mysimdata,rho.list,c=c,d=d)
      }else{mcmcobj3 <- mcmcobj2}
      mcmcobj4 <- Update_Mixture_NSBSR(mcmcobj3,mysimdata)
      
      new.mcmcobj <- mcmcobj4
      
      if(iter>n.burn){
        
        new.mcmcobj$EstiMat[(iter-n.burn),] <- c(new.mcmcobj$Parset$beta,
                                                 new.mcmcobj$Parset$tau.e,
                                                 new.mcmcobj$Parset$theta,
                                                 new.mcmcobj$Parset$rho.th1,
                                                 new.mcmcobj$Parset$rho.th2)
        if((iter-n.burn)>1){
          
          Esti.summary = summary(mcmc(new.mcmcobj$EstiMat[1:(iter-n.burn),]))
          
          PredResult <- Predict_NSBSR(new.mcmcobj,Esti.summary,mysimdata)
          
          # source("Suppl/Show_MCMC_NSBSR_Images.R")
          
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
  
  save.image(paste0("MODIS/",Sys.Date(),"_NSBSR(nd=",nd,")).RData"))
  
}



Y = Yl[1:29,1:29]

Pred.result = summary(mcmc.list(NSBSR.PredMat))
n = nrow(Pred.result$statistics)/3
temp = matrix(Pred.result$statistics[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.result$statistics[(n+1):(2*n),1],ncol=sqrt(n))
Yhat = temp[1:29,1:29] + temp2[1:29,1:29]

par(mfrow=c(1,2))
fields::image.plot(Y)
fields::image.plot(Yhat)

sdhat = sqrt(Pred.result$statistics[1:n,2]^2 + Pred.result$statistics[(1+n):(2*n),2]^2)
sdhatmat = matrix(sdhat,sqrt(n),sqrt(n))
mean((Y > Yhat - 1.96 * sdhatmat[1:29,1:29]) * (Y < Yhat + 1.96 * sdhatmat[1:29,1:29]))


