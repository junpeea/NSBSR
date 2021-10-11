
rm(list=ls())
setwd("~/Desktop/Code3(Acting)")
load("200601StudyDataset.RData")
# load("2018AURAdataset.RData")
# load("2017AURAdataset.RData")

X = as.data.frame(X.pred); colnames(X) <- c("1","R","S","NO2","VOC","UV","SO2");head(X)
Y = matrix(mydata$O3L3,ncol=1); colnames(Y) <- "O3L3"
simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
               X=X, Y=Y)
mysimdata = data.frame(Y,X,simdata$map)

par(mfrow=c(1,2))
fields::image.plot(matrix(simdata$Y,(xrange+1),(yrange+1)))
# fields::image.plot(matrix(simdata$Y,(xrange+1),(yrange+1))[seq(1,xrange+1,by=2),seq(1,yrange+1,by=2)])
fields::image.plot(matrix(as.matrix(simdata$X) %*% fit1$coefficients,(xrange+1),(yrange+1)))

#######################################################################################################################
# Lattice Kriging
#######################################################################################################################

library(LatticeKrig)
setwd("~/Desktop/heatoncomparison/Code/LatticeKrig")
source("LKrig.sim.conditional.foreach.R")
library(spam64)
library(dplyr)

time.LK <- system.time({
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  # train.sim_data$z = fit0$residuals
  # train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  formula <- O3L3~R+S+NO2+VOC+UV+SO2
  
  O3.loc = cbind(train.sim_data$x,train.sim_data$y)
  O3.concentration = train.sim_data$O3L3
  O3.covariate = train.sim_data %>% dplyr::select(R,S,NO2,VOC,UV,SO2) %>% cbind(1) %>% as.matrix
  obj  <- LatticeKrig(O3.loc, O3.concentration, z=O3.covariate)
  
  M <- 250
  outputSim<- LKrig.sim.conditional.foreach(obj, nx=(xrange+1), ny=(yrange+1), M = M, nCores=4)
  
  standardError<- sqrt(apply( outputSim$g.draw, 1, "var") + obj$sigma.MLE^2)
  int = 2*1.96*standardError
  
  pred.Lkrig <- predict(obj,xnew=cbind(mysimdata$x,mysimdata$y),znew=X.pred)
  pred.Lkrig.mat <- matrix(pred.Lkrig,(xrange+1),(yrange+1))
  
})

par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.Lkrig.mat,zlim=c(270,350))

mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.Lkrig.mat))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.Lkrig.mat))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.Lkrig.mat) )^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.Lkrig.mat) )^2,na.rm=T));round(rmse,3)
round(mean(int),3)
time.LK

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(standardError,(xrange+1),(yrange+1))
mean((Realdata > pred.Lkrig.mat - 1.96 * sdhatmat) * (Realdata < pred.Lkrig.mat + 1.96 * sdhatmat))

par(mfrow=c(1,2))
fields::image.plot(Realdata)
fields::image.plot(pred.Periodic)

# setwd("~/Desktop/Code3(Acting)")
# save.image(paste0(Sys.Date(),"_LKrig.Rdata"))
# save.image(paste0(Sys.Date(),"_LKrig2018.Rdata"))
# save.image(paste0(Sys.Date(),"_LKrig2017.Rdata"))

# pdf(paste0(Sys.Date(),"_AURA_LK.pdf"),width=7*1.5,height=7*1.5)
# pdf(paste0(Sys.Date(),"_AURA_LK2018.pdf"),width=7*1.5,height=7*1.5)
# pdf(paste0(Sys.Date(),"_AURA_LK2017.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(pred.Lkrig.mat,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat",zlim=c(270,350))
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
# dev.off()

