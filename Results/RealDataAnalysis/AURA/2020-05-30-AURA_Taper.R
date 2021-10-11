
rm(list=ls())
setwd("~/Desktop/Code3(Acting)")
load("200601StudyDataset.RData")

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
# Tapering
#######################################################################################################################

library(dplyr)
library(spam)
library(fields)

time.taper02 <- system.time({
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  Esti.grid = train.sim_data %>% dplyr::select("x","y");dim(Esti.grid)
  Pred.grid = mydata  %>% dplyr::select("x","y");dim(Pred.grid)
  # Xobs <- train.sim_data %>% dplyr::select(R,S,NO2,VOC,UV,SOl,SO2) %>% cbind(1) %>% as.matrix
  Xobs <- train.sim_data %>% dplyr::select(R,S,NO2,VOC,UV,SO2) %>% cbind(1) %>% as.matrix
  Yobs <- train.sim_data %>% dplyr::select(O3L3) %>% as.matrix
  
  temp = Pred.grid[,c("x","y")]
  Rstar = matrix(c(1,0,0,1),2,2)
  an.loc  = t(Rstar%*%t(temp))
  distvec1 = as.vector(as.matrix(dist(an.loc)))
  temp = Esti.grid[,c("x","y")]
  Rstar = matrix(c(1,0,0,1),2,2)
  an.loc  = t(Rstar%*%t(temp))
  distvec2 = as.vector(as.matrix(dist(an.loc)))
  
  delta = 0.2*max(Esti.grid)
  hobs <- nearest.dist(Esti.grid, miles=FALSE, delta=delta)
  hobs <- hobs + t(hobs)
  C_taper = cov.wend1( hobs, theta=c( delta,1,0))
  
  predcormtx = matrix(Matern(distvec1,smoothness=0.5,range=10),ncol=nrow(Pred.grid))
  esticormtx = matrix(Matern(distvec2,smoothness=0.5,range=10),ncol=nrow(Esti.grid))
  
  C_tapered = esticormtx * C_taper
  inv_cov_n = solve(C_tapered)
  H.mtx = predcormtx[,!is.na(mydata$O3L3)]
  
  M = 250
  # esticovmtx = summary(fit0)$sigma*esticormtx
  esticovmtx = esticormtx
  cholSig = chol(esticovmtx)
  esticor = as.spam(x=esticovmtx)
  sample = rmvnorm.spam(n=M,Sigma = esticor,Rstruct=cholSig)
  sample.Y = apply(sample,1,function(w){
    as.matrix(X)%*%matrix(fit0$coefficients,ncol=1) + H.mtx %*% inv_cov_n %*% w  
  } )
  
  pred.taper20 = as.matrix(X)%*%matrix(fit0$coefficients,ncol=1) + H.mtx %*% inv_cov_n %*% fit0$residuals
  int = (apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.975)
         -apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.025))
  pred.taper20 = matrix(pred.taper20,(xrange+1),(yrange+1))
  
})

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(standardError,(xrange+1),(yrange+1))
mean((Realdata > apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.025)) * (Realdata < apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.975)))


par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.taper20,zlim=c(270,350))

mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.taper20))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.taper20))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.taper20))^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.taper20))^2,na.rm=T));round(rmse,3)
round(mean(int),3)
time.taper02

setwd("~/Desktop/Code3(Acting)")
save.image(paste0(Sys.Date(),"_Taper.Rdata"))

pdf(paste0(Sys.Date(),"_AURA_Taper.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(pred.taper20,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()

