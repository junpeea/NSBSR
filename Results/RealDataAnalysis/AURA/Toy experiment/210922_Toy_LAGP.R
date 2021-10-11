
#######################################################################################################################
# LAGP
#######################################################################################################################

library(LatticeKrig)
library(tgp)
library(laGP)
library(plgp)

time.LAGP <- system.time({
  
  n <- 30
  temp = expand.grid(c(0:(n-1)),c(0:(n-1)))
  sim_process = data.frame(x=temp$Var1,y=temp$Var2)
  sim_process$Y = as.vector(Y0)
  nr <- 30; nc <- 30
  Lon <- matrix(sim_process$x, nrow=nr, ncol=nc)
  Lat <- matrix(sim_process$y, nrow=nr, ncol=nc)
  
  ind = is.na(as.vector(Y00));summary(ind)
  ## extract the training/testing problem
  nas <- which(ind)
  X <- cbind(sim_process$x, sim_process$y)[-nas,] # training input
  N <- nrow(X)
  y <- sim_process$Y[-nas] # training resposne
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
  
  ind = is.na(as.vector(Y00));summary(ind)
  sim_data <- sim_process[which(ind==FALSE),]
  sim_data$x1 = rep(1,nrow(sim_data))
  sim_data$x2 = scale(sim_process$x[which(ind==FALSE)])
  sim_data$x3 = scale(sim_process$y[which(ind==FALSE)])
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  
  pred.LAGP = as.vector(Y00)
  pred.LAGP[ind] <- psub$mean
  pred.LAGP = matrix(pred.LAGP,30,30)
  
  sdhatvec = rep(0,900)
  sdhatvec[ind] = stdev
  sdhatmat = matrix(sdhatvec,30,30)
  mean((Realdata > pred.LAGP - 1.96 * sdhatmat) * (Realdata < pred.LAGP + 1.96 * sdhatmat))+0.25
  
  
  fields::image.plot(pred.LAGP)
  
  pred.NSBSR = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  
  summary((t(mylist[[k]])[1:30,1:30] - t(pred.LAGP[1:30,1:30]) )^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.LAGP[1:30,1:30])%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.LAGP[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
  summary((t(mylist[[k]])[1:30,1:30] - t(pred.NSBSR[1:30,1:30]))^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.NSBSR[1:30,1:30])%>%as.vector,use="complete.obs")
  # summary((t(mylist[[1]])[1:30,1:30] - t(temp3[1:30,1:30]) )^2 %>% as.vector)
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.NSBSR[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(temp3[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
})

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.LAGP[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.LAGP[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int),4)
time.LAGP

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_LAGP(toy).Rdata"))

