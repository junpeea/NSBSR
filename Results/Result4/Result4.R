########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Script name : NSBSR_Result4

# Date : 201908xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################

############################################################################################## SETTINGS
rm(list = ls());
# setwd("C:/Users/user/Desktop/180102_LCY_THESIS") # set your main directory
setwd("D:/NSBSR_FINAL/Library")
source(file="PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="NSBSR_Func.R")             # Functions for Nonarametric Spectral Bayesian Spatial Regression
######################################## Available versions
# source(file="NSBSR_Func_MCMC_ver1.R")      # Version1. Isotropic theta
# source(file="NSBSR_Func_MCMC_ver2.R")        # Version2. common version
# source(file="NSBSR_Func_MCMC_ver2_Rcpp.R")   # Version2. common version
# source(file="NSBSR_Func_MCMC_ver3.R")      # Version3. 2-sided rho-hyper-par
########################################
require(rgdal)
require(gdalUtils)
require(raster)
require(ggmap)


# ### LOAD MODIS dataset ##############################################################################################################
# setwd("C:/Users/user/Desktop/180102_LCY_THESIS/MODIS_PROJECT")
# sds <- list()
# gdalinfo("MOD08_E3.A2017145.006.2017153180952.hdf"); sds[[1]] <- get_subdatasets("MOD08_E3.A2017145.006.2017153180952.hdf")
# gdalinfo("MOD08_E3.A2017153.006.2017163121114.hdf"); sds[[2]] <- get_subdatasets("MOD08_E3.A2017153.006.2017163121114.hdf")
# gdalinfo("MOD08_E3.A2017161.006.2017171120513.hdf"); sds[[3]] <- get_subdatasets("MOD08_E3.A2017161.006.2017171120513.hdf")
# gdalinfo("MOD08_E3.A2017169.006.2017178120811.hdf"); sds[[4]] <- get_subdatasets("MOD08_E3.A2017169.006.2017178120811.hdf")
# gdalinfo("MOD08_E3.A2017177.006.2017192120526.hdf"); sds[[5]] <- get_subdatasets("MOD08_E3.A2017177.006.2017192120526.hdf")
# gdalinfo("MOD08_E3.A2017185.006.2017193082915.hdf"); sds[[6]] <- get_subdatasets("MOD08_E3.A2017185.006.2017193082915.hdf")
# gdalinfo("MOD08_E3.A2017193.006.2017201083533.hdf"); sds[[7]] <- get_subdatasets("MOD08_E3.A2017193.006.2017201083533.hdf")
# gdalinfo("MOD08_E3.A2017201.006.2017212170536.hdf"); sds[[8]] <- get_subdatasets("MOD08_E3.A2017201.006.2017212170536.hdf")
# gdalinfo("MOD08_E3.A2017209.006.2017217085735.hdf"); sds[[9]] <- get_subdatasets("MOD08_E3.A2017209.006.2017217085735.hdf")
# gdalinfo("MOD08_E3.A2017217.006.2017234125534.hdf"); sds[[10]] <- get_subdatasets("MOD08_E3.A2017217.006.2017234125534.hdf")
# gdalinfo("MOD08_E3.A2017225.006.2017233084637.hdf"); sds[[11]] <- get_subdatasets("MOD08_E3.A2017225.006.2017233084637.hdf")
# gdalinfo("MOD08_E3.A2017233.006.2017249151735.hdf"); sds[[12]] <- get_subdatasets("MOD08_E3.A2017233.006.2017249151735.hdf")
# gdalinfo("MOD08_E3.A2017241.006.2017249084348.hdf"); sds[[13]] <- get_subdatasets("MOD08_E3.A2017241.006.2017249084348.hdf")
# #####################################################################################################################################

### LOAD MODIS dataset ##############################################################################################################
setwd("C:/Users/user/Desktop/NSBSR_FINAL/MODIS")
sds <- list()
gdalinfo("MOD08_D3.A2019203.061.2019204194413.hdf"); sds[[1]] <- get_subdatasets("MOD08_D3.A2019203.061.2019204194413.hdf")
gdalinfo("MOD08_D3.A2019204.061.2019205090003.hdf"); sds[[2]] <- get_subdatasets("MOD08_D3.A2019204.061.2019205090003.hdf")
gdalinfo("MOD08_D3.A2019205.061.2019206085102.hdf"); sds[[3]] <- get_subdatasets("MOD08_D3.A2019205.061.2019206085102.hdf")
gdalinfo("MOD08_D3.A2019206.061.2019207090537.hdf"); sds[[4]] <- get_subdatasets("MOD08_D3.A2019206.061.2019207090537.hdf")
gdalinfo("MOD08_D3.A2019207.061.2019208091228.hdf"); sds[[5]] <- get_subdatasets("MOD08_D3.A2019207.061.2019208091228.hdf")
gdalinfo("MOD08_D3.A2019208.061.2019209092216.hdf"); sds[[6]] <- get_subdatasets("MOD08_D3.A2019208.061.2019209092216.hdf")
gdalinfo("MOD08_D3.A2019209.061.2019210090322.hdf"); sds[[7]] <- get_subdatasets("MOD08_D3.A2019209.061.2019210090322.hdf")
gdalinfo("MOD08_D3.A2019210.061.2019211085540.hdf"); sds[[8]] <- get_subdatasets("MOD08_D3.A2019210.061.2019211085540.hdf")
gdalinfo("MOD08_D3.A2019211.061.2019212090039.hdf"); sds[[9]] <- get_subdatasets("MOD08_D3.A2019211.061.2019212090039.hdf")
#####################################################################################################################################


mylist = list()
colindex = rainbow(77,alpha=0.5)
for(dayind in 1:9){
  
  gdal_translate(sds[[dayind]][844], dst_dataset = paste0("Ozone_","QAmean_",dayind,".tif"))
  rast <- raster(paste0("Ozone_","QAmean_",dayind,".tif"))
  totdata <- as.matrix(rast);dim(totdata)
  mylist[[dayind]] <- log(totdata[(24+90):(53+90),(112+180):(141+180)]+1e-10)
  
}

# image plot over googlemap
register_google(key="AIzaSyCjtzlHwGhbI0qqwF5AukJ_LcjfZ0P-ffY")
mymap  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "satellite");plot(mymap)
mymap2  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "terrain-labels");plot(mymap2)
mymap3 = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "toner-lines");plot(mymap3)


for(dayind in 1:9){
  pdf(paste0(Sys.Date(),"MODDay",dayind,".pdf"))
  par(mfrow=c(1,1),mar=c(1,1,1,5))
  x = seq(112,141,by=1);y = 28/28*(x-113)+25
  plot(x,y,type="n",xlab="lon",ylab="lat",main="",xaxt="n",yaxt="n")
  axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
  axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
  fields::image.plot(x=x,y=seq(24,53,by=1),z=t(mylist[[dayind]]),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
  par(new=T);plot(mymap2,xaxt="n",xaxt="n",yaxt="n")
  par(new=T);plot(mymap3,xaxt="n",xaxt="n",yaxt="n")
  dev.off()
}

#####################################################################################################################################
setwd("C:/Users/user/Desktop/NSBSR_FINAL")
source("Library/NSBSR_ImportLibrary.R")
source("Library/NSBSR_Func.R")

# NSBSR VALIDATION
start = Sys.time()

for(nd in 3:3){
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  
  Realdata = mylist[[nd]]
  
  # Real data (example)
  par(mfrow=c(1,1))
  x = seq(112,141,by=1);y = 28/28*(x-113)+25
  plot(x,y,type="n",xlab="lon",ylab="lat",main="")
  image(x=x,y=seq(24,53,by=1),z=t(Realdata),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
  par(new=T);plot(mymap2)
  par(new=T);plot(mymap3)
  # dev.off()
  
  grid.size = c(length(x),length(y)); obs.delta = c(2,2)               # Environment settings
  ND <- 1
  X.seed <- 1; Y.seed <- 1  
  
  sp1 = seq(0,(grid.size[1]-1),by=1); sp2 = seq(0,(grid.size[2]-1),by=1)
  space.grid = expand.grid(sp1,sp2)
  omega.grid = expand.grid(2*pi*sp1/grid.size[1],2*pi*sp2/grid.size[2])
  
  nloc = nrow(space.grid)
  
  bddind = rep(1,nloc)
  bddind[which(space.grid[,1] %% (grid.size[1]/2) == 0 )] <- 2
  bddind[which(space.grid[,2] %% (grid.size[2]/2) == 0 )] <- 2
  plot(space.grid,col=bddind,main="Show_Boundaries")
  
  obsind = matrix(NA,nrow=nrow(t(Realdata)),ncol=ncol(t(Realdata)))
  obsind[seq(1,nrow(obsind),by=2),seq(1,ncol(obsind),by=2)] <- 4
  obsind[is.na(t(Realdata))%>%as.vector & obsind==4] <- 2
  obsind = obsind %>% as.vector
  nobs = length(which(obsind > 0))
  obs.prob = length(which(obsind==4))/nobs
  plot(space.grid,col=obsind,main="Show_MonitoredLocation")
  
  map.grid = data.frame(space.grid,omega.grid,bddind,obsind); colnames(map.grid) <- c("x","y","w.x","w.y","bdd","obs")
  
  d.o.o.w   <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.x","w.y")) %>% dist %>% as.matrix
  d.o.o.w.x <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.x")) %>% dist %>% as.matrix
  d.o.o.w.y <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.y")) %>% dist %>% as.matrix
  
  # source("Common/Get2dGridcovariates.R")                # Generate 2-dimension lattice grid covariates
  X = cbind(1,expand.grid((x-mean(x))/sd(x),(y-mean(y))/sd(y))); dim(X); colnames(X) <- c("X0","X1","X2")
  Y = matrix(t(Realdata),ncol=1); dim(Y)
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 X=X, Y=Y)
  
  par(mfrow=c(1,1))
  fields::image.plot(matrix(simdata$Y,30,30),col=colindex[6:68])
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=1000; d=1
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  # source("Common/Initialize_NSBSR.R")
  source("Common/Initialize_INBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=1000,n.burn=9)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  # source("Library/NSBSR_Func_MCMC_ver7.R")
  source("Library/NSBSR_Func_MCMC_ver9.R")
  # source("Library/NSBSR_Func_MCMC_ver10.R")
  # source("Library/NSBSR_Func_MCMC_ver11.R")
  
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
          
          # source("Suppl/Show_MCMC_NSBSR_Images.R")
          
          par(mfrow=c(1,3))
          X = simdata$X; Y = simdata$Y
          
          Yhat = PredResult$result1 + PredResult$result2
          brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
          score = sqrt(mean((Y-Yhat)^2))
          
          Ymat  = matrix(Y,   ncol=sqrt(simdata$nloc))
          Yhmat = matrix(Yhat,ncol=sqrt(simdata$nloc))
          
          # Take 1
          fields::image.plot(Ymat,main="TRUE",col=colindex[6:68])
          # Take 2
          fields::image.plot(Yhmat,main=paste0("MSE=",round(score,2)),col=colindex[6:68])
          
          plot(Ymat,Yhmat);abline(b=1,lty=2,col=2)
          
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
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result4")
  save.image(paste0(Sys.Date(),"_MODISresult1_nd=",nd,".Rdata"))
  
}

record = Sys.time() - start
record


# PBSR validation

start = Sys.time()
for(nd in 3:3){
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  
  Realdata = mylist[[nd]]
  
  grid.size = c(length(x),length(y)); obs.delta = c(2,2)               # Environment settings
  X.seed <- 1; Y.seed <- 1  
  
  sp1 = seq(0,(grid.size[1]-1),by=1); sp2 = seq(0,(grid.size[2]-1),by=1)
  space.grid = expand.grid(sp1,sp2)
  omega.grid = expand.grid(2*pi*sp1/grid.size[1],2*pi*sp2/grid.size[2])
  
  nloc = nrow(space.grid)
  
  bddind = rep(1,nloc)
  bddind[which(space.grid[,1] %% (grid.size[1]/2) == 0 )] <- 2
  bddind[which(space.grid[,2] %% (grid.size[2]/2) == 0 )] <- 2
  plot(space.grid,col=bddind,main="Show_Boundaries")
  
  obsind = matrix(NA,nrow=nrow(t(Realdata)),ncol=ncol(t(Realdata)))
  obsind[seq(1,nrow(obsind),by=2),seq(1,ncol(obsind),by=2)] <- 4
  obsind[is.na(t(Realdata))%>%as.vector & obsind==4] <- 2
  obsind = obsind %>% as.vector
  nobs = length(which(obsind > 0))
  obs.prob = length(which(obsind==4))/nobs
  plot(space.grid,col=obsind,main="Show_MonitoredLocation")
  
  map.grid = data.frame(space.grid,omega.grid,bddind,obsind); colnames(map.grid) <- c("x","y","w.x","w.y","bdd","obs")
  X = cbind(1,expand.grid((x-mean(x))/sd(x),(y-mean(y))/sd(y))); dim(X); colnames(X) <- c("X0","X1","X2")
  Y = matrix(t(Realdata),ncol=1); dim(Y)
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 X=X, Y=Y)
  
  par(mfrow=c(1,1))
  fields::image.plot(matrix(simdata$Y,30,30),col=colindex[6:68])
  
  ################################################################################## Initialize PBSR ###
  source("Common/Initialize_PBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=100,n.burn=9)
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
    n.iter = PBSR00.mcmcobj0$n.iter
    n.burn = PBSR00.mcmcobj0$n.burn
    
    iter <- 1
    for(iter in 1:n.iter){
      
      PBSR00.new.mcmcobj <- Gibbs_PBSR(PBSR00.mcmcobj0,simdata,fix.kappa=FALSE,a=10,b=50)
     
      if(iter>n.burn){
        
        PBSR00.new.mcmcobj$EstiMat[(iter-n.burn),] <- c(PBSR00.new.mcmcobj$Parset$beta,
                                                        PBSR00.new.mcmcobj$Parset$kappa,
                                                        PBSR00.new.mcmcobj$Parset$nu2,
                                                        PBSR00.new.mcmcobj$Parset$sigma2,
                                                        PBSR00.new.mcmcobj$Parset$range)

        PredResult <- Predict_PBSR(PBSR00.new.mcmcobj,simdata)
        PBSR00.new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                                        as.vector(PredResult$result2),
                                                        as.vector(PredResult$result3))
        
        # # PredResult <- Predict_PBSR(PBSR01.new.mcmcobj,simdata)
        # #######################################################################
        # # setwd("C:/Users/user/Desktop/WORK2018/180527_NSBSR")
        # # png(paste0("temp_chain",chain,"iter",iter,".png"))
        par(mfrow=c(1,3))
        X = simdata$X; Y = simdata$Y
        Yhat = PredResult$result1 + PredResult$result2
        brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
        score = sqrt(mean((Y-Yhat)^2))
        
        Ymat =  matrix(Y,   ncol=sqrt(simdata$nloc))
        Yhmat = matrix(Yhat,ncol=sqrt(length(Yhat)))
        
        col <- rev(rainbow(77)[1:63])
        
        fields::image.plot(Ymat,main="TRUE",col=colindex[6:68])
        fields::image.plot(Yhmat,main=paste0("MSE=",round(score,2)),col=colindex[6:68])
        # fields::image.plot(Yhat,breaks=brk,col=col,main=round(score,2))
        plot(Ymat,Yhmat);abline(b=1,lty=2,col=2)
        # par(mfrow=c(1,1))
        # fields::image.plot(PredResult$result4)
        
        # # dev.off()
 
        
      }
      
      PBSR00.mcmcobj0 <- PBSR00.new.mcmcobj
      
      cat("PBSR_chain",chain,"/",n.chain,"iter",iter,"/",n.iter,
          "alpha",PBSR00.new.mcmcobj$Parset$kappa," ",date(),'\n')
      
    }
    
    PBSR00.EstiMat[[chain]] <- PBSR00.new.mcmcobj$EstiMat
    PBSR00.PredMat[[chain]] <- PBSR00.new.mcmcobj$PredMat
    
  }
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result4")
  save.image(paste0(Sys.Date(),"_MODISresult2_nd=",nd,".Rdata"))
  
}

record = Sys.time() - start
record


mse_r2_table = matrix(NA,9,4)
for(k in 1:9){
  # setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result4")
  setwd("D:/NSBSR_FINAL/THESIS_Results(New)/Result4")
  load(file=paste0("2019-08-06_MODISresult1_nd=",k,".Rdata"))
  Pred.summary = summary(mcmc.list(NSBSR.PredMat))
  n = nrow(Pred.summary$statistics)/3
  temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
  temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
  temp3 = temp + temp2
  
  sc = sqrt(mean(t(mylist[[k]])[1:30,1:30] - t(temp3[1:30,1:30]),na.rm=T)^2)
  r2 = summary(lm(as.vector(t(mylist[[k]])[1:30,1:30]) ~ as.vector(t(temp3))))$r.squared
  mse_r2_table[k,1]<-round(sc,4)
  mse_r2_table[k,3]<-round(r2,2)
  
  load(file=paste0("2019-08-06_MODISresult2_nd=",k,".Rdata"))
  Pred.summary = summary(mcmc.list(PBSR00.PredMat))
  
  n = nrow(Pred.summary$statistics)/3
  temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
  temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
  temp3 = temp + temp2
  
  sc = sqrt(mean(t(mylist[[k]])[1:30,1:30] - t(temp3[1:30,1:30]),na.rm=T)^2)
  r2 = summary(lm(as.vector(t(mylist[[k]])[1:30,1:30]) ~ as.vector(t(temp3))))$r.squared
  mse_r2_table[k,2]<-round(sc,4)
  mse_r2_table[k,4]<-round(r2,2)
}
mse_r2_table

plot(mse_r2_table[,1],col=4,type='o',ylim=c(0,0.20))
par(new=T);plot(mse_r2_table[,2],col=2,type='o',ylim=c(0,0.20))






# NSBSR APPLICATION
start = Sys.time()
for(nd in 3:3){
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  
  Realdata = mylist[[nd]]
  
  # Real data
  colindex = rainbow(77,alpha=0.5)
  par(mfrow=c(1,1))
  x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
  plot(x,y,type="n",xlab="lon",ylab="lat",main="")
  Z = matrix(NA,59,59)
  Z[seq(1,59,by=2),seq(1,59,by=2)] = t(Realdata)
  # Z[seq(2,59,by=2),seq(2,59,by=2)] = t(Realdata)
  image(x=x,y=seq(24,53,by=0.5),z=Z,col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
  par(new=T);plot(mymap2)
  par(new=T);plot(mymap3)
  
  setwd("C:/Users/user/Desktop/NSBSR_FINAL")
  grid.size = c(length(x),length(y)); obs.delta = c(2,2)               # Environment settings
  X.seed <- 1; Y.seed <- 1  
  
  sp1 = seq(0,(grid.size[1]-1),by=1); sp2 = seq(0,(grid.size[2]-1),by=1)
  space.grid = expand.grid(sp1,sp2)
  omega.grid = expand.grid(2*pi*sp1/grid.size[1],2*pi*sp2/grid.size[2])
  
  nloc = nrow(space.grid)
  
  bddind = rep(1,nloc)
  bddind[which(space.grid[,1] %% (grid.size[1]/2) == 0 )] <- 2
  bddind[which(space.grid[,2] %% (grid.size[2]/2) == 0 )] <- 2
  plot(space.grid,col=bddind,main="Show_Boundaries")
  
  obsind = rep(NA,nloc)
  obsind[which(rowSums(space.grid %% obs.delta) == 0 )] <- 2
  obsind[which(is.na(Z)==0)] <- 4
  nobs = length(which(is.na(obsind)==0))
  # set.seed(1); p = rbinom(nobs,1,obs.prob)
  plot(space.grid,col=obsind,main="Show_MonitoredLocation")
  
  map.grid = data.frame(space.grid,omega.grid,bddind,obsind); colnames(map.grid) <- c("x","y","w.x","w.y","bdd","obs")
  
  d.o.o.w   <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.x","w.y")) %>% dist %>% as.matrix
  d.o.o.w.x <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.x")) %>% dist %>% as.matrix
  d.o.o.w.y <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.y")) %>% dist %>% as.matrix
  
  
  X = cbind(1,expand.grid((x-mean(x))/sd(x),(y-mean(y))/sd(y))); dim(X); colnames(X) <- c("X0","X1","X2")
  Y = matrix(Z,ncol=1); dim(Y)
  simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
                 X=X, Y=Y)
  # table(simdata$map$obs,useNA = c("always"))
  
  par(mfrow=c(1,1))
  fields::image.plot(matrix(simdata$Y,59,59),col=colindex[6:68])
  
  # Set NSBSR Hyperparameters
  a=100; b=100; c=1000; d=100
  rho.vec  = exp(seq(log(1e-05),log(1e-03),length.out=10))
  rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
  rho.list = list(rho.vec,rho.grid)
  
  ################################################################################## Initialize NSBSR ###
  # source("Common/Initialize_NSBSR.R")
  source("Common/Initialize_INBSR.R")
  null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=100,n.burn=9)
  n.chain = length(null.mcmcobj.list)
  NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
  init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                    function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )
  
  ############################################ MCMC - NSBSR #################################################
  # source("Library/NSBSR_Func_MCMC_ver7.R")
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
          
          par(mfrow=c(1,3),mar=c(2,2,3,3))
          X = simdata$X; Y = simdata$Y
          
          Yhat = PredResult$result1 + PredResult$result2
          brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
          score = sqrt(mean((Y-Yhat)^2))
          
          Ymat  = matrix(Y,   ncol=sqrt(simdata$nloc))
          Yhmat = matrix(Yhat,ncol=sqrt(simdata$nloc))
          
          col <- rev(rainbow(77)[1:63])
          
          # Take 1
          fields::image.plot(Ymat,main="TRUE",axes=F)
          # Take 2
          fields::image.plot(Yhmat,main=paste0("MSE=",round(score,2)),axes=F)
          # Take 3
          plot(Y,Yhat);abline(b=1,lty=2,col=2,axes=F)
          
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
  
  # setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result4")
  # save.image(paste0(Sys.Date(),"MODISresult3.Rdata"))
  
}

record = Sys.time() - start
record







k <- 3
pdf(paste0(Sys.Date(),"MODIS_FIG1.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
brk = (seq(0.99*min(mylist[[k]],na.rm=T),1.01*max(mylist[[k]],na.rm=T),length.out=78))
x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[k]])),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()

pdf(paste0(Sys.Date(),"MODIS_FIG2.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
brk = (seq(0.99*min(mylist[[k]],na.rm=T),1.01*max(mylist[[k]],na.rm=T),length.out=78))
x = seq(112,141,by=2);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="TRAINING DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=2),z=(t(mylist[[k]][seq(1,30,by=2),seq(1,30,by=2)])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()


Pred.summary = summary(mcmc.list(NSBSR.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2

k <- 3
pdf(paste0(Sys.Date(),"MODIS_FIG3.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
brk = (seq(0.99*min(mylist[[k]],na.rm=T),1.01*max(mylist[[k]],na.rm=T),length.out=78))
x = seq(112,141,by=1);y = 28/28*(x-113)+25
sc = mean(t(mylist[[3]])[1:30,1:30] - t(temp3),na.rm=T)^2
r2 = summary(lm(as.vector(t(mylist[[3]])[1:30,1:30]) ~ as.vector(t(temp3))))$r.squared
plot(x,y,type="n",xlab="",ylab="",main=paste0("PREDICTED (NSBSR) OZONE"),xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,54,by=5),labels=seq(25,54,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()



Pred.summary = summary(mcmc.list(PBSR00.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2


pdf(paste0(Sys.Date(),"MODIS_FIG4.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
brk = (seq(0.99*min(mylist[[3]],na.rm=T),1.01*max(mylist[[3]],na.rm=T),length.out=78))
x = seq(112,141,by=1);y = 28/28*(x-113)+25
sc = mean((mylist[[3]] - t(temp3))^2,na.rm=T)
r2 = summary(lm(as.vector(mylist[[3]]) ~ as.vector(t(temp3))))$r.squared
par(mfrow=c(1,1),mar=c(2,2,2,5))
plot(x,y,type="n",xlab="",ylab="",main=paste0("Predicted (PBSR) OZONE"),xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(24,53,by=5),labels=seq(24,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()




getwd()

colindex = rainbow(77,alpha=0.5)
brk = (seq(7.8,8.4,length.out=78))

# x = seq(112,141,by=1);y = 28/28*(x-113)+25
# plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n')
# title(xlab="lon",ylab="lat",line=2)
# axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
# axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
# fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[1]])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
# par(new=T);plot(mymap2)
# par(new=T);plot(mymap3)


Pred.summary = summary(mcmc.list(NSBSR.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,3],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),3],ncol=sqrt(n))
temp3 = temp + temp2

summary(temp3 %>% as.vector)

pdf(paste0(Sys.Date(),"MODIS_FIG5.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=temp3,col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()



n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2

pdf(paste0(Sys.Date(),"MODIS_FIG6.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE (2.5% quantile)",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=temp3,col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,5],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),5],ncol=sqrt(n))
temp3 = temp + temp2


pdf(paste0(Sys.Date(),"MODIS_FIG7.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE (97.5% quantile)",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=temp3,col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()

