


#######################################################################################################################
# FRK
#######################################################################################################################

library(sp)
library(FRK)

time.FRK <- system.time({
  
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
  coordinates(sim_data) = ~x + y # change into an sp object
  
  ## Make BAUs as SpatialPixels
  
  grid_BAUs <- sim_process                   # assign BAUs
  
  # BAUs$Missing <- is.na(BAUs$Temp)      # mark which BAUs contain missing data
  
  grid_BAUs$Y <- NULL                     # remove data from BAUs
  
  grid_BAUs$fs <- 1                          # set fs variation to unity
  
  coordinates(grid_BAUs)  <- ~x + y        # convert to SpatialPointsDataFrame
  
  gridded(grid_BAUs) <- TRUE                 # convert to SpatialPixelsDataFrame
  
  
  # grid_BAUs <- auto_BAUs(manifold=plane(),data=sim_data,
  #                        nonconvex_hull=FALSE,cellsize = c(1.33),type="grid")
  # grid_BAUs$fs = 1
  
  G <- auto_basis(manifold = plane(),
                  data=sim_data,
                  nres = 3,
                  regular = 1,
                  type = "bisquare",
                  subsamp = 20000)
  
  f <- z ~ 1
  
  S <- FRK(f = f,                       # formula for SRE model
           
           data = sim_data,             # data
           
           basis = G,               # Basis
           
           BAUs = grid_BAUs,                 # BAUs
           
           tol = 0.01)                   # EM iterations
  
  
  ### Check fit info
  
  
  ### Predict over BAUs
  grid_BAUs <- predict(S)
  
  ## Predict
  # BAUs_pred <- SRE.predict(S)           # predict over all BAUs
  BAUs_pred <- predict(S)           # predict over all BAUs
  BAUs_pred_df <- data.frame(BAUs_pred) # convert to data frame
  
  ## Compute variance of predicted observations and conf. intervals
  BAUs_pred_df$sd_obs <- sqrt(BAUs_pred_df$sd^2 + S@Ve[1,1])
  int <- 2*1.96*BAUs_pred_df$sd_obs
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  pred.FRK   = matrix(X%*%fit0$coefficients + grid_BAUs$mu,30,30)

  par(mfrow=c(1,2))
  fields::image.plot(Realdata)
  fields::image.plot(pred.FRK)
  
  sdhat = BAUs_pred_df$sd_obs
  sdhatmat = matrix(sdhat,n,n)
  mean((Realdata > pred.FRK - 1.96 * sdhatmat) * (Realdata < pred.FRK + 1.96 * sdhatmat))
  
  
  # pred.NSBSR = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  # 
  # summary((t(mylist[[K]])[1:30,1:30] - t(pred.FRK[1:30,1:30]) )^2 %>% as.vector)
  # cor(x=t(mylist[[K]])[1:30,1:30]%>%as.vector,y=t(pred.FRK[1:30,1:30])%>%as.vector,use="complete.obs")
  # par(mfrow=c(2,2))
  # fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(pred.FRK[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # 
  # summary((t(mylist[[K]])[1:30,1:30] - t(pred.NSBSR[1:30,1:30]))^2 %>% as.vector)
  # cor(x=t(mylist[[K]])[1:30,1:30]%>%as.vector,y=t(pred.NSBSR[1:30,1:30])%>%as.vector,use="complete.obs")
  # # summary((t(mylist[[1]])[1:30,1:30] - t(temp3[1:30,1:30]) )^2 %>% as.vector)
  # par(mfrow=c(2,2))
  # fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(pred.NSBSR[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # # fields::image.plot(temp3[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
})

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.FRK[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.FRK[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int),4)
time.FRK

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_FRK(toy).Rdata"))


