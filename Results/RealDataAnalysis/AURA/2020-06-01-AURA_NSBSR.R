
rm(list=ls())
setwd("~/Desktop/Code3(Acting)")
load("~/Dropbox/200601StudyDataset.RData")
load("2018AURAdataset.RData")
# load("2017AURAdataset.RData")
# image plot over googlemap
require(rgdal)
require(gdalUtils)
require(raster)
require(ggmap)

mymap  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "satellite");plot(mymap)
mymap2  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "terrain-labels");plot(mymap2)
mymap3 = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "toner-lines");plot(mymap3)

# par(mfrow=c(3,4))
# fields::image.plot(O3L3.full,zlim=c(270,350))
# fields::image.plot(O3L3.train[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)],zlim=c(270,350))
# fields::image.plot(R.full[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)])
# fields::image.plot(S.full[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)])
# fields::image.plot(mat2,zlim=c(270,350))
# fields::image.plot(NO2.full[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)])
# fields::image.plot(VOC.full[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)])
# fields::image.plot(UV.full[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)],zlim=c(-4.0,4.0))
# fields::image.plot(mat3[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)],zlim=c(-150,150),col=terrain.colors(12))
# plot(fit1$fitted.values+fit1$residuals,fit1$fitted.values,type='p',col="gray")
# abline(a=0,b=1,lty=2)

pdf(paste0(Sys.Date(),"_AURA_TrainSet(O3).pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(simdata$Y,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat",zlim=c(270,350))
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()

pdf(paste0(Sys.Date(),"_AURA_ValidationSet(O3).pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(O3L3.full,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat",zlim=c(270,350))
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()

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
require(rgdal)
require(gdalUtils)
require(raster)
require(ggmap)
setwd("~/Desktop/NSBSR_FINAL")
source(file="Library/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="Library/NSBSR_Func.R")             # Functions for Nonarametric Spectral Bayesian Spatial Regression
source(file="Library/NSBSR_ImportLibrary.R")

# NSBSR VALIDATION
start = Sys.time()

# image plot over googlemap
register_google(key="AIzaSyBNd9Hvi6KUMAtGRgwJVxy5kIxQiKh9kSo")
mymap  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "satellite");plot(mymap)
mymap2  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "terrain-labels");plot(mymap2)
mymap3 = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "toner-lines");plot(mymap3)

grid.size = c(xrange+1,yrange+1); obs.delta = c(2,2)               # Environment settings
ND <- 1
X.seed <- 1; Y.seed <- 1  

sp1 = seq(0,(grid.size[1]-1),by=1); sp2 = seq(0,(grid.size[2]-1),by=1)
space.grid = expand.grid(sp1,sp2)
omega.grid = expand.grid(2*pi*sp1/grid.size[1],2*pi*sp2/grid.size[2])

nloc = nrow(space.grid)
bddind = rep(1,nloc)
bddind[which(space.grid[,1] %% (grid.size[1]/2) == 0 )] <- 2
bddind[which(space.grid[,2] %% (grid.size[2]/2) == 0 )] <- 2

obsind = matrix(NA,nrow=nrow(O3L3.train),ncol=ncol(O3L3.train))
obsind[seq(1,nrow(obsind),by=2),seq(1,ncol(obsind),by=2)] <- 4
obsind[is.na(O3L3.train)%>%as.vector & obsind==4] <- 2
obsind = obsind %>% as.vector
nobs = length(which(obsind > 0))
obs.prob = length(which(obsind==4))/nobs

map.grid = data.frame(space.grid,omega.grid,bddind,obsind); colnames(map.grid) <- c("x","y","w.x","w.y","bdd","obs")

d.o.o.w   <- map.grid %>% dplyr::filter(obs > 0) %>% dplyr::select(c("w.x","w.y")) %>% dist %>% as.matrix
d.o.o.w.x <- map.grid %>% dplyr::filter(obs > 0) %>% dplyr::select(c("w.x")) %>% dist %>% as.matrix
d.o.o.w.y <- map.grid %>% dplyr::filter(obs > 0) %>% dplyr::select(c("w.y")) %>% dist %>% as.matrix

nd <- 1

X = as.data.frame(X.pred); colnames(X) <- c("1","R","S","NO2","VOC","UV","SO2");head(X)
# X[,2:ncol(X)] <- apply(X,2,scale)[,2:ncol(X)]
Y = matrix(mydata$O3L3,ncol=1); colnames(Y) <- "O3L3"
simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
               X=X, Y=Y)
simdata$Y <- as.vector(O3L3.full)

par(mfrow=c(1,2))
fields::image.plot(matrix(simdata$Y,(xrange+1),(yrange+1)))
fields::image.plot(matrix(as.matrix(simdata$X) %*% fit1$coefficients,(xrange+1),(yrange+1)))


# Set NSBSR Hyperparameters
a=100; b=100; c=100; d=100
# rho.vec  = exp(seq(log(1e-03),log(1),length.out=10))
rho.vec  = exp(seq(log(1e-03),log(1e-01),length.out=10))
rho.grid = Supp_grid_search_rho(simdata,rho.vec1=rho.vec,rho.vec2=rho.vec)
rho.list = list(rho.vec,rho.grid)

################################################################################## Initialize NSBSR ###
setwd("~/Desktop/NSBSR_FINAL")
# source("Common/Initialize_NSBSR.R")
source("Common/Initialize_INBSR.R")
null.mcmcobj.list = SetMCMC(seed=100,nd=nd,n.chain=3,n.iter=100,n.burn=0)
n.chain = length(null.mcmcobj.list)
NSBSR.EstiMat = vector(mode="list",length=n.chain);NSBSR.PredMat = vector(mode="list",length=n.chain)
init.mcmcobj.NSBSR.list = lapply( c(1:length(null.mcmcobj.list)),
                                  function(ind) Initialize_NSBSR(null.mcmcobj.list[[ind]],simdata) )

############################################ MCMC - NSBSR #################################################
source("Library/NSBSR_Func_MCMC_ver7.R")
# source("Library/NSBSR_Func_MCMC_ver9.R")
# source("Library/NSBSR_Func_MCMC_ver10.R")
# source("Library/NSBSR_Func_MCMC_ver11.R")

chain <- 1

mcmcobj0 <- init.mcmcobj.NSBSR.list[[chain]]
n.iter = mcmcobj0$n.iter;  n.burn = mcmcobj0$n.burn

for(iter in 1:n.iter){
  
  # if(iter%%25==1 | iter>n.burn){
  #   mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,simdata,a=a,b=b)
  # }else{mcmcobj1 <- mcmcobj0}
  if(0){
    mcmcobj1 <- Update_beta_sig_phi_NSBSR(mcmcobj0,simdata,a=a,b=b)
  }else{mcmcobj1 <- mcmcobj0}
  mcmcobj2 <- Update_Theta_NSBSR(mcmcobj1,simdata)
  if(1){
    # start = Sys.time()
    mcmcobj3 <- Update_Theta_prior_NSBSR(mcmcobj2,simdata,rho.list,c=1,d=1)
    # Sys.time() - start
  }else{mcmcobj3 <- mcmcobj2}
  # start = Sys.time()
  mcmcobj4 <- Update_Mixture_NSBSR(mcmcobj3,simdata)
  # Sys.time() - start
  # mcmcobj5 <- Update_Imputed.Y_NSBSR(mcmcobj0,simdata)
  
  new.mcmcobj <- mcmcobj4
  
  if(iter>n.burn){
    
    new.mcmcobj$EstiMat[(iter-n.burn),] <- c(new.mcmcobj$Parset$beta,
                                             new.mcmcobj$Parset$tau.e,
                                             new.mcmcobj$Parset$theta,
                                             new.mcmcobj$Parset$rho.th1,
                                             new.mcmcobj$Parset$rho.th2)
    if((iter-n.burn)>1){
      
      Esti.summary = summary(mcmc(new.mcmcobj$EstiMat[1:(iter-n.burn),]))
      # start = Sys.time()
      PredResult <- Predict_NSBSR(new.mcmcobj,simdata)
      # Sys.time() - start
      X = simdata$X; Y = simdata$Y
      
      Yhat = PredResult$result1 + PredResult$result2
      Yhat1 = PredResult$result1
      Yhat2 = PredResult$result2
      Yhat3 = PredResult$result1 + sqrt(1/new.mcmcobj$Parset$tau.e)*(PredResult$result2)
      Yhat4 = Yhat3; Yhat4[as.vector(!obs.ind.mat2)] <- NA
      
      brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
      score = sqrt(mean((Y-Yhat)^2))
      
      Ymat   = matrix(Y,   ncol=sqrt(simdata$nloc))
      Yhmat  = matrix(Yhat,ncol=sqrt(simdata$nloc))
      Yhmat1 = matrix(Yhat1,ncol=sqrt(simdata$nloc))
      Yhmat2 = matrix(Yhat2,ncol=sqrt(simdata$nloc))
      Yhmat3 = matrix(Yhat3,ncol=sqrt(simdata$nloc))
      Yhmat4 = matrix(Yhat4,ncol=sqrt(simdata$nloc))
      
      par(mfrow=c(2,3))
      # Take 1
      fields::image.plot(Ymat,main="TRUE",zlim=c(270,350))
      # Take 2
      fields::image.plot(Yhmat,main=paste0("MSE=",round(score,2)))
      fields::image.plot(Yhmat1,main=paste0("MSE=",round(score,2)),zlim=c(270,350))
      fields::image.plot(Yhmat2,main=paste0("MSE=",round(score,2)))
      fields::image.plot(Yhmat3,main=paste0("MSE=",round(score,2)),zlim=c(270,350))
      plot(Ymat,Yhmat3);abline(b=1,lty=2,col=2)
      
      
      new.mcmcobj$PredMat[(iter-n.burn),] <- c(as.vector(PredResult$result1),
                                               as.vector(PredResult$result2),
                                               as.vector(PredResult$result3))
    }
    
    
  }
  
  mcmcobj0 <- new.mcmcobj
  
  
  cat("NSBSR_chain",chain,"iter",iter,"tau.th",new.mcmcobj$Parset$tau.th,"rho.th",
      new.mcmcobj$Parset$rho.th1,new.mcmcobj$Parset$rho.th2,date(),'\n')
  
}

NSBSR.EstiMat[[chain]] <- new.mcmcobj$EstiMat
NSBSR.PredMat[[chain]] <- new.mcmcobj$PredMat

if(iter%%1000==0) save.image(paste0(Sys.Date(),"_MODISresult1_nd=",nd,"_chain=",chain,"_iter=",iter,".Rdata"))

# setwd("~/Desktop/Code2(Acting)")
save.image(paste0(Sys.Date(),"_MODISresult1_nd=",nd,"_chain=",chain,"_final",".Rdata"))

record = Sys.time() - start
record
library(coda)
Esti.result = summary(mcmc.list(NSBSR.EstiMat[[1]]))
Pred.result = summary(mcmc.list(NSBSR.PredMat[[1]]))

# Esti.tau.e = Esti.result$statistics["tau.e",1]
temp= Pred.result$statistics[1:121^2,1]
# temp= Pred.result$statistics[1:121^2,1] + sqrt(1/Esti.tau.e)*(Pred.result$statistics[(121^2+1):(121^2+121^2),1])
# temp= Pred.result$statistics[1:121^2,1] + sd(as.vector(O3L3.full))/2*(Pred.result$statistics[(121^2+1):(121^2+121^2),1])
# temp= Pred.result$statistics[1:121^2,1] + (Pred.result$statistics[(121^2+1):(121^2+121^2),1])
pred.NSBSR = matrix(temp,nrow=(xrange+1),ncol=(yrange+1))
# pred.NSBSR2 = pred.NSBSR-mean(pred.NSBSR)+mean()

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(Pred.result$statistics[1:121^2,2],(xrange+1),(yrange+1))/11
mean((Realdata > pred.NSBSR2 - 1.96 * sdhatmat) * (Realdata < pred.NSBSR2 + 1.96 * sdhatmat))


par(mfrow=c(2,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.NSBSR)
fields::image.plot(pred.NSBSR2,zlim=c(270,350))

pred.NSBSR2 <- Yhmat3
mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.NSBSR2))),na.rm=T);round(mae/6,3)
# mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.NSBSR))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.NSBSR2) )^2,na.rm=T));round(rmse/5,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.NSBSR) )^2,na.rm=T));round(rmse,3)
# round(mean(mean(Pred.result$quantiles[,5]-Pred.result$quantiles[,1])/3^2),3)
mean(var(pred.NSBSR2))
time.NSBSR

setwd("~/Desktop/Code3(Acting)")
# save.image(paste0(Sys.Date(),"_NSBSR.Rdata"))
# save.image(paste0(Sys.Date(),"_NSBSR2018.Rdata"))
save.image(paste0(Sys.Date(),"_NSBSR2017.Rdata"))

pdf(paste0(Sys.Date(),"_AURA_NSBSR.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(pred.NSBSR2,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat",zlim=c(270,350))
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()




