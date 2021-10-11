
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
# Spartial Partitioning
#######################################################################################################################

library(fields)
library(dplyr)

time.partition4 <- system.time({
  
  sim_process <- mysimdata
  sim_process$area = rep(NA,nrow(sim_process))
  ind = list()
  ind[[1]] = which(sim_process$x<mean(sim_process$x)  & sim_process$y<mean(sim_process$y))
  ind[[2]] = which(sim_process$x<mean(sim_process$x)  & sim_process$y>=mean(sim_process$y))
  ind[[3]] = which(sim_process$x>=mean(sim_process$x) & sim_process$y<mean(sim_process$y))
  ind[[4]] = which(sim_process$x>=mean(sim_process$x) & sim_process$y>=mean(sim_process$y))
  sim_process$area[ind[[1]]] <- 1
  sim_process$area[ind[[2]]] <- 2
  sim_process$area[ind[[3]]] <- 3
  sim_process$area[ind[[4]]] <- 4
  sim_process$obs = !is.na(mydata$O3L3)
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  Esti.grid = train.sim_data   %>% dplyr::select("x","y");dim(Esti.grid)
  Pred.grid = sim_process  %>% dplyr::select("x","y");dim(Pred.grid)
  # Xobs <- train.sim_data %>% dplyr::select(R,S,NO2,VOC,UV,SOl,SO2) %>% cbind(1) %>% as.matrix
  Xobs <- train.sim_data %>% dplyr::select(R,S,NO2,VOC,UV,SO2) %>% cbind(1) %>% as.matrix
  Yobs <- train.sim_data %>% dplyr::select(O3L3) %>% as.matrix
  
  X = X.pred
  
  fit = list(); Esti.grid = list(); Pred.grid = list()
  pred.partition    = rep(NA,nrow(sim_process))
  M = 250; sample.Y = matrix(NA,nrow(sim_process),ncol=M)
  # i <- 1
  for(i in 1:4){
    sim_process_i = (sim_process %>% filter(area==i));dim(sim_process_i)
    train.ind_i = sim_process_i$obs
    sim_data_i = sim_process_i[train.ind_i,]
    # fit[[i]] = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=sim_data_i)
    fit[[i]] = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=sim_data_i)
    Esti.grid[[i]] = sim_data_i  %>% dplyr::select("x","y"); dim(Esti.grid[[i]])
    Pred.grid[[i]] = sim_process_i                %>% dplyr::select("x","y"); dim(Pred.grid[[i]])
    
    temp = Pred.grid[[i]][,c("x","y")]
    distvec1 = as.vector(as.matrix(dist(temp)))
    temp = Esti.grid[[i]][,c("x","y")]
    distvec2 = as.vector(as.matrix(dist(temp)))
    
    predcormtx = matrix(Matern(distvec1,smoothness=0.5,range=10),ncol=nrow(Pred.grid[[i]]))
    esticormtx = matrix(Matern(distvec2,smoothness=0.5,range=10),ncol=nrow(Esti.grid[[i]]))
    inv_cov_n = solve(esticormtx);dim(inv_cov_n)
    H.mtx = predcormtx[,train.ind_i]; dim(H.mtx)
    
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
  
  
  
  pred.partition4 = matrix(pred.partition,(xrange+1),(yrange+1))
  int = (apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.975)
         -apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.025))
})

par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.partition4 ,zlim=c(270,350))

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(standardError,(xrange+1),(yrange+1))
mean((Realdata >= apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.025)) * (Realdata <= apply(sample.Y[which(rowSums(is.na(sample.Y))==0),],1,quantile, probs=0.975)))


mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.partition4))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.partition4))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.partition4) )^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.partition4) )^2,na.rm=T));round(rmse,3)
round(mean(int),3)
time.partition4


setwd("~/Desktop/Code3(Acting)")
save.image(paste0(Sys.Date(),"_Partition.Rdata"))

