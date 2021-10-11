#######################################################################################################################
# Bayesian Nearest Neighbor Gausian Process
#######################################################################################################################

# rm(list=ls())
library(spNNGP)
library(MBA)
library(fields)

time.nngp <- system.time({
  
  n <- 30
  temp = expand.grid(c(0:(n-1)),c(0:(n-1)))
  sim_process = data.frame(x=temp$Var1,y=temp$Var2)
  sim_process$Y = as.vector(Y0)
  
  scale = function(x) {(x-mean(x))/sd(x)}
  
  ind = is.na(as.vector(Y00));summary(ind)
  sim_data <- sim_process[which(ind==FALSE),]
  sim_data$x1 = rep(1,nrow(sim_data))
  sim_data$x2 = scale(sim_process$x[which(ind==FALSE)])
  sim_data$x3 = scale(sim_process$y[which(ind==FALSE)])
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  sim_data$z = fit0$residuals
  sim_data$std <- 0.1
  head(sim_data)
  
  
  coords = cbind(sim_data$x,sim_data$y)
  cov.model <- "exponential"
  sigma.sq <- 6.5
  sigma.sq.IG <- c(2, sigma.sq)
  
  
  g <- 5
  theta.alpha <- as.matrix(expand.grid(seq(7, 9, length.out=g), seq(0.00001/sigma.sq, 0.001/sigma.sq, length.out=g)))
  
  colnames(theta.alpha) <- c("phi", "alpha")
  
  
  m.c <- spConjNNGP(Y~x1+x2+x3-1, data=sim_data, coords=coords, n.neighbors = 15,
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
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X.ho = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  coords.ho = cbind(sim_process$x,sim_process$y)
  
  m.p <- spConjNNGP(Y~x1+x2+x3-1, data=sim_data, coords=coords, n.neighbors = 15,
                    X.0=X.ho, coords.0=coords.ho,
                    n.omp.threads = 10,
                    theta.alpha = theta.alpha, sigma.sq.IG = sigma.sq.IG,
                    cov.model = cov.model)
  
  m.p$run.time[3]/60
  
  pred.nngp.mat = matrix(m.p$y.0.hat,30,30)
  
  a <- m.p$ab[2]
  t <- qt(0.975, 2*a)
  
  me <- sqrt((a-1)/a*m.p$y.0.hat.var)*t
  
  pred <- cbind(m.p$y.0.hat, m.p$y.0.hat-me, m.p$y.0.hat+me)
  
  colnames(pred) <- c("50%","2.5%","97.5%")
  
  sdhatmat = matrix(standardError,n,n)
  mean((Realdata > matrix(pred[,2],30,30)) * (Realdata < matrix(pred[,3],30,30)))
  
  
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  pred.NSBSR = matrix(X.ho%*%fit0$coefficients + as.vector(temp2),30,30)
  
  
  summary((t(mylist[[k]])[1:30,1:30] - t(pred.nngp.mat[1:30,1:30]) )^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.nngp.mat[1:30,1:30])%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.nngp.mat[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
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

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.nngp.mat[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.nngp.mat[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
int  = 2*me; round(mean(int),4)
time.nngp

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_NNGP(toy).Rdata"))



