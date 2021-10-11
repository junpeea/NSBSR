#######################################################################################################################
# Covariance Tapering
#######################################################################################################################

time.taper02 <- system.time({
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
  head(sim_process)
  Esti.grid = sim_data     %>% dplyr::select("x","y");dim(Esti.grid)
  Pred.grid = sim_process  %>% dplyr::select("x","y");dim(Pred.grid)
  Xobs <- sim_data %>% dplyr::select(x1,x2,x3) %>% as.matrix
  Yobs <- sim_data %>% dplyr::select(Y) %>% as.matrix
  
  temp = Pred.grid[,c("x","y")]
  Rstar = matrix(c(1,0,0,1),2,2)
  an.loc  = t(Rstar%*%t(temp))
  distvec1 = as.vector(as.matrix(dist(an.loc)))
  temp = Esti.grid[,c("x","y")]
  Rstar = matrix(c(1,0,0,1),2,2)
  an.loc  = t(Rstar%*%t(temp))
  distvec2 = as.vector(as.matrix(dist(an.loc)))
  library(spam)
  delta = 0.2*max(Esti.grid)
  hobs <- nearest.dist(Esti.grid, miles=FALSE, delta=delta)
  hobs <- hobs + t(hobs)
  C_taper = cov.wend1( hobs, theta=c( delta,1,0))
  
  library(fields)
  predcormtx = matrix(Matern(distvec1,smoothness=0.5,range=10),ncol=nrow(Pred.grid))
  esticormtx = matrix(Matern(distvec2,smoothness=0.5,range=10),ncol=nrow(Esti.grid))
  
  C_tapered = esticormtx * C_taper
  inv_cov_n = solve(C_tapered)
  H.mtx = predcormtx[,!ind]
  
  M = 250
  # esticovmtx = summary(fit0)$sigma*esticormtx
  esticovmtx = esticormtx
  cholSig = chol(esticovmtx)
  esticor = as.spam(x=esticovmtx)
  sample = rmvnorm.spam(n=M,Sigma = esticor,Rstruct=cholSig)
  sample.Y = apply(sample,1,function(w){
    X%*%fit0$coefficients + H.mtx %*% inv_cov_n %*% w  
  } )
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  pred.taper20 = X%*%fit0$coefficients + H.mtx %*% inv_cov_n %*% fit0$residuals
  int = (apply(sample.Y,1,quantile, probs=0.975) - apply(sample.Y,1,quantile, probs=0.025))/sqrt(M)
  pred.taper20 = matrix(pred.taper20,30,30)
  
  sdhatmat = matrix(standardError,n,n)
  mean((Realdata > apply(sample.Y,1,quantile, probs=0.025)) * (Realdata < apply(sample.Y,1,quantile, probs=0.975)))
  
  pred.NSBSR = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  
  summary((t(mylist[[k]])[1:30,1:30] - t(pred.taper20[1:30,1:30]) )^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.taper20[1:30,1:30])%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.taper20[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
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

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.taper20[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.taper20[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int),4)
time.taper02

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_Taper(toy).Rdata"))

