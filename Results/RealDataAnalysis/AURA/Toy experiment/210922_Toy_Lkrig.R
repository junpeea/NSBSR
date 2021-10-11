#######################################################################################################################
# Lattice Kriging
#######################################################################################################################

library(LatticeKrig)
setwd("C:/Users/yjun/Documents/GitHub/heatoncomparison/Code/LatticeKrig")
source("LKrig.sim.conditional.foreach.R")

time.LK <- system.time({
  ### Generate process and data
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
  # coordinates(sim_data) = ~x + y # change into an sp object
  
  O3.loc = cbind(sim_data$x,sim_data$y)
  O3.concentration = sim_data$Y
  O3.covariate = cbind(sim_data$x2,sim_data$x3)
  # O3.covariate = sim_data$x3
  obj  <- LatticeKrig(O3.loc, O3.concentration,z=O3.covariate)
  
  M <- 250
  outputSim<- LKrig.sim.conditional.foreach(obj, nx=30, ny=30, M = M, nCores=2)
  
  names(outputSim)
  
  standardError<- sqrt(apply( outputSim$g.draw, 1, "var") + obj$sigma.MLE^2)
  int = 2*1.96*standardError
  
  # compare with and without linear covariates
  pred.Lkrig <- predict(obj,xnew=cbind(sim_process$x,sim_process$y))
  pred.Lkrig.mat <- matrix(pred.Lkrig,30,30)
  
  sdhatmat = matrix(standardError,n,n)
  mean((Realdata > pred.FRK - 1.96 * sdhatmat) * (Realdata < pred.FRK + 1.96 * sdhatmat))
  
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  pred.NSBSR = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  
  summary((t(mylist[[k]])[1:30,1:30] - t(pred.Lkrig.mat)[1:30,1:30])^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.Lkrig.mat)%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.Lkrig.mat[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
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

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.Lkrig.mat[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.Lkrig.mat[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int),4)
time.LK


setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_LK(toy).Rdata"))


