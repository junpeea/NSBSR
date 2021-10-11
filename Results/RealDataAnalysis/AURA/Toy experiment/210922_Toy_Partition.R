#######################################################################################################################
# Spartial Partitioning
#######################################################################################################################
library(fields)

time.partition4 <- system.time({
  
  ### Generate process and data
  n <- 30
  temp = expand.grid(c(0:(n-1)),c(0:(n-1)))
  sim_process = data.frame(x=temp$Var1,y=temp$Var2)
  sim_process$Y = as.vector(Y0)
  scale = function(x) {(x-mean(x))/sd(x)}
  sim_process$area = rep(NA,nrow(sim_process))
  ind = list()
  ind[[1]] = which(sim_process$x<15 & sim_process$y<15)
  ind[[2]] = which(sim_process$x<15 & sim_process$y>=15)
  ind[[3]] = which(sim_process$x>=15 & sim_process$y<15)
  ind[[4]] = which(sim_process$x>=15 & sim_process$y>=15)
  ind0 = is.na(as.vector(Y00));summary(ind0)
  sim_process$area[ind[[1]]] <- 1
  sim_process$area[ind[[2]]] <- 2
  sim_process$area[ind[[3]]] <- 3
  sim_process$area[ind[[4]]] <- 4
  obsind = !is.na(as.vector(Y00));summary(obsind)
  sim_process$obs = as.numeric(obsind)
  
  sim_data <- sim_process[obsind,];dim(sim_data)
  sim_data$x1 = rep(1,nrow(sim_data))
  sim_data$x2 = scale(sim_process$x[obsind])
  sim_data$x3 = scale(sim_process$y[obsind])
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  # Xobs <- sim_data %>% dplyr::select(x1,x2,x3) %>% as.matrix
  # Yobs <- sim_data %>% dplyr::select(Y) %>% as.matrix
  
  fit = list(4); Esti.grid = list(4); Pred.grid = list(4)
  mx  = mean(sim_process$x[which(ind0==FALSE)]) ;my  = mean(sim_process$y[which(ind0==FALSE)])
  sdx = sd(sim_process$x[which(ind0==FALSE)])   ;sdy = sd(sim_process$y[which(ind0==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  
  pred.partition    = rep(NA,nrow(sim_process))
  M = 250; sample.Y = matrix(NA,nrow(sim_process),ncol=M)
  for(i in 1:4){
    sim_process_i = (sim_process %>% filter(area==i));dim(sim_process_i)
    ind_i = sim_process_i$obs==1;summary(ind_i)
    sim_data_i    = (sim_data %>% filter(area==i));dim(sim_data_i)
    fit[[i]] = lm(Y~-1+x1+x2+x3,data=sim_data_i)
    Esti.grid[[i]] = sim_data_i          %>% dplyr::select("x","y")
    Pred.grid[[i]] = sim_process_i       %>% dplyr::select("x","y")
    
    temp = Pred.grid[[i]][,c("x","y")]
    distvec1 = as.vector(as.matrix(dist(temp)))
    temp = Esti.grid[[i]][,c("x","y")]
    distvec2 = as.vector(as.matrix(dist(temp)))
    
    predcormtx = matrix(Matern(distvec1,smoothness=0.5,range=10),ncol=nrow(Pred.grid[[i]]))
    esticormtx = matrix(Matern(distvec2,smoothness=0.5,range=10),ncol=nrow(Esti.grid[[i]]))
    inv_cov_n = solve(esticormtx);dim(inv_cov_n)
    H.mtx = predcormtx[,ind_i]; dim(H.mtx)
    
    pred.partition[ind[[i]]] <- X[ind[[i]],]%*%fit[[i]]$coefficients + H.mtx %*% inv_cov_n %*% fit[[i]]$residuals
    
    esticovmtx = esticormtx; dim(esticovmtx)
    cholSig = chol(esticovmtx)
    esticov = as.spam(x=esticovmtx)
    sample = rmvnorm.spam(n=M,Sigma = esticov,Rstruct=cholSig);dim(sample)
    sample.Y_i = apply(sample,1,function(w){
      X[ind[[i]],]%*%fit[[i]]$coefficients + H.mtx %*% inv_cov_n %*% w  
    } )
    sample.Y[ind[[i]],] <- sample.Y_i
  }
  
  pred.partition4 = matrix(pred.partition,30,30)
  pred.NSBSR = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  int = (apply(sample.Y,1,quantile, probs=0.975) - apply(sample.Y,1,quantile, probs=0.025))/sqrt(M)
  
  mean((Realdata > apply(sample.Y,1,quantile, probs=0.025)) * (Realdata < apply(sample.Y,1,quantile, probs=0.975)))
  
  summary((t(mylist[[1]])[1:30,1:30] - t(pred.partition4[1:30,1:30]) )^2 %>% as.vector)
  cor(x=t(mylist[[1]])[1:30,1:30]%>%as.vector,y=t(pred.partition4[1:30,1:30])%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.partition4[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
  summary((t(mylist[[1]])[1:30,1:30] - t(pred.NSBSR[1:30,1:30]))^2 %>% as.vector)
  cor(x=t(mylist[[1]])[1:30,1:30]%>%as.vector,y=t(pred.NSBSR[1:30,1:30])%>%as.vector,use="complete.obs")
  # summary((t(mylist[[1]])[1:30,1:30] - t(temp3[1:30,1:30]) )^2 %>% as.vector)
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.NSBSR[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(temp3[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
})

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.partition4[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.partition4[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int),4)
time.partition4


setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_Partition(toy).Rdata"))

