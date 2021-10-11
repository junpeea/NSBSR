
rm(list=ls())
setwd("~/Desktop/Code3(Acting)")
# load("200601StudyDataset.RData")
# load("2018AURAdataset.RData")
load("2017AURAdataset.RData")

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
# LAGP
#######################################################################################################################

library(LatticeKrig)
library(tgp)
library(laGP)
library(plgp)

time.LAGP <- system.time({
  
  sim_process <- mysimdata
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  nr <- nx; nc <- ny
  Lon <- matrix(sim_process$x, nrow=nr, ncol=nc)
  Lat <- matrix(sim_process$y, nrow=nr, ncol=nc)
  
  ind = is.na(mydata$O3L3);summary(ind)
  ## extract the training/testing problem
  nas <- which(ind)
  X <- cbind(sim_process$x, sim_process$y)[-nas,] # training input
  N <- nrow(X)
  y <- sim_process$O3L3[-nas] # training resposne
  XX <- cbind(sim_process$x, sim_process$y)[nas,] # test input
  NN <- nrow(XX)
  
  ## code inputs and predictive grid to unit cube
  maxX <- apply(rbind(X, XX), 2, max)
  minX <- apply(rbind(X, XX), 2, min)
  for (j in 1:ncol(X)){
    X[,j] <- X[,j] - minX[j]	
    X[,j] <- X[,j]/(maxX[j]-minX[j])  
    XX[,j] <- XX[,j] - minX[j]	
    XX[,j] <- XX[,j]/(maxX[j]-minX[j])
  }
  
  ##
  ## time for fitting
  ##
  
  ## macro-scale analysis on a maximum entropy sub-design of size n=100
  n <- 100
  sub <- dopt.gp(100, Xcand=X)
  
  ## priors for the global (subset) GP
  da <- darg(list(mle=TRUE, max=10), X)
  ga <- garg(list(mle=TRUE, max=10), y)
  
  
  ## fit the global GP
  gpsepi <- newGPsep(sub$XX, y[sub$fi], d=da$start, g=ga$start, dK=TRUE)
  that <- mleGPsep(gpsepi, param="both", tmin=c(da$min, ga$min), 
                   tmax=c(da$max, ga$max), ab=c(da$ab, ga$ab), maxit=200)
  
  ## predictions from the global GP on the test set
  psub <- predGPsep(gpsepi, XX, lite=TRUE)
  
  ## calculation of residuals at the full set of input locations
  pX <- predGPsep(gpsepi, X, lite=TRUE)
  yresid <- y - pX$mean
  
  ## Now for laGP on residuals from the (subset) global GP
  ## scale the inputs according to the macro-analysis lengthscales
  scale <- sqrt(that$theta[1:2])
  Xs <- X; XXs <- XX
  for(j in 1:ncol(Xs)){
    Xs[,j] <- Xs[,j] / scale[j]
    XXs[,j] <- XXs[,j] / scale[j]
  }
  
  ## local analysis on residuals (and scaled inputs) from local analysis
  out <- aGPsep(Xs, yresid, XXs, d=list(start=1, max=20), g=that$theta[3], 
                omp.threads=4, verb=0)
  
  ## Deriving the uncertainty in the global GP mean requires 
  ## removing the nugget from psub variance.
  K <- covar.sep(sub$XX, d=that$theta[1:2], g=that$theta[3])
  psi <- drop(t(y[sub$fi]) %*% solve(K) %*% y[sub$fi])
  s2.f <- psub$s2 * nrow(sub$XX) / psi
  eps <- sqrt(.Machine$double.eps)
  s2.f <- s2.f - that$theta[3] + eps
  s2.f <- s2.f * psi / nrow(sub$XX)
  
  ## combine variances from psub (from residuals) and local GP
  stdev <- sqrt(s2.f + out$var)
  int <- 2*1.96*stdev
  
  pred.LAGP = vector(len=nrow(sim_process))
  pred.LAGP[-nas]  <- train.sim_data$O3L3
  pred.LAGP[nas] <- psub$mean
  pred.LAGP = matrix(pred.LAGP,(xrange+1),(yrange+1))
  
})

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(0,(xrange+1),(yrange+1))
sdhatmat[ind] <- stdev 
mean((Realdata >= pred.LAGP - 1.96 * sdhatmat) * (Realdata <= pred.LAGP + 1.96 * sdhatmat))


par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.LAGP,zlim=c(270,350))

mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.LAGP))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.LAGP))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.LAGP) )^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.LAGP) )^2,na.rm=T));round(rmse,3)
round(mean(int),3)
time.LAGP

setwd("~/Desktop/Code3(Acting)")
# save.image(paste0(Sys.Date(),"_LAGP.Rdata"))
# save.image(paste0(Sys.Date(),"_LAGP2018.Rdata"))
save.image(paste0(Sys.Date(),"_LAGP2017.Rdata"))

# pdf(paste0(Sys.Date(),"_AURA_LAGP.pdf"),width=7*1.5,height=7*1.5)
pdf(paste0(Sys.Date(),"_AURA_LAGP2018.pdf"),width=7*1.5,height=7*1.5)
# pdf(paste0(Sys.Date(),"_AURA_LAGP2017.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(pred.LAGP,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat",zlim=c(270,350))
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()
