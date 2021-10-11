
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
# Bayesian Nearest Neighbor Gausian Process
#######################################################################################################################

library(spNNGP)
library(MBA)
library(fields)

time.nngp <- system.time({
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  coords = cbind(train.sim_data$x,train.sim_data$y)
  cov.model <- "exponential"
  sigma.sq <- 6.5
  sigma.sq.IG <- c(2, sigma.sq)
  
  
  g <- 5
  theta.alpha <- as.matrix(expand.grid(seq(7, 9, length.out=g), seq(0.00001/sigma.sq, 0.001/sigma.sq, length.out=g)))
  colnames(theta.alpha) <- c("phi", "alpha")
  
  # formula <- O3L3~R+S+NO2+VOC+UV+SOl+SO2
  formula <- O3L3~R+S+NO2+VOC+UV+SO2
  m.c <- spConjNNGP(formula, data=train.sim_data, coords=coords, n.neighbors = 15,
                    k.fold = 5, score.rule = "crps",
                    n.omp.threads = 3,
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                    cov.model = cov.model)
  
  
  m.c$run.time[3]/60
  
  m.c$sigma.sq.hat
  tau.sq <- m.c$theta.alpha[2]*m.c$sigma.sq.hat
  tau.sq
  
  ##prediction
  theta.alpha <- as.vector(m.c$theta.alpha)
  names(theta.alpha) <- c("phi", "alpha")
  
  X.ho <- X.pred
  coords.ho = cbind(mydata$x,mydata$y)
  
  m.p <- spConjNNGP(formula, data=train.sim_data, coords=coords, n.neighbors = 15,
                    X.0=X.ho, coords.0=coords.ho,
                    n.omp.threads = 10,
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                    cov.model = cov.model)
  
  m.p$run.time[3]/60
  
  pred.nngp.mat = matrix(m.p$y.0.hat,(xrange+1),(yrange+1))
  
  a <- m.p$ab[2]
  t <- qt(0.975, 2*a)
  
  me <- sqrt((a-1)/a*m.p$y.0.hat.var)*t
  
  pred <- cbind(m.p$y.0.hat, m.p$y.0.hat-me, m.p$y.0.hat+me)
  
  colnames(pred) <- c("50%","2.5%","97.5%")
  
})

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(me,(xrange+1),(yrange+1))
mean((Realdata > pred.nngp.mat - sdhatmat) * (Realdata < pred.nngp.mat + sdhatmat))


par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.nngp.mat,zlim=c(270,350))

mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.nngp.mat))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.nngp.mat))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.nngp.mat) )^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.nngp.mat) )^2,na.rm=T));round(rmse,3)
int  = 2*me; round(mean(int,na.rm=T),3)
time.nngp

setwd("~/Desktop/Code3(Acting)")
save.image(paste0(Sys.Date(),"_NNGP.Rdata"))

pdf(paste0(Sys.Date(),"_AURA_NNGP.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(pred.nngp.mat,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat",zlim=c(270,350))
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()


