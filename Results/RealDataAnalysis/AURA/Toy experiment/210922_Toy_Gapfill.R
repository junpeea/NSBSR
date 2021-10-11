#######################################################################################################################
# Gapfill
#######################################################################################################################

## all packages are available on CRAN
library("gapfill") 
library("doParallel");
registerDoParallel(2) # run 40 tasks in parallel
library("abind")

time.Gapfill <- system.time({
  n <- 30
  temp = expand.grid(c(0:(n-1)),c(0:(n-1)))
  sim_process = data.frame(x=temp$Var1,y=temp$Var2)
  sim_process$Y = as.vector(Y0)
  
  scale = function(x) {(x-mean(x))/sd(x)}
  
  ind = is.na(as.vector(Y00));summary(ind)
  # sim_data <- sim_process[which(ind==FALSE),]
  sim_data <- sim_process
  sim_data$x1 = rep(1,nrow(sim_data))
  sim_data$x2 = scale(sim_process$x)
  sim_data$x3 = scale(sim_process$y)
  sim_data$Y  = as.vector(Y00)
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  
  ## rearrange data as matrix
  data <- with(sim_data,
               array(Y, c(length(unique(x)), length(unique(y)))))
  dim <- dim(data)
  nx <- dim[1]; ny <- dim[2]
  
  ## display data
  ## Image(data)
  
  ## augment data: since the gapfill method is designed for
  ## spatio-temproal data, we artificially create 9 additional
  ## and similar images by shifting the given image
  data_augmented <- abind(data,
                          data[c(1,1:(nx-1)),],
                          data[,c(1,1:(ny-1))],
                          data[c(2:nx,nx),],
                          data[,c(2:ny,ny)],
                          
                          data[c(2:nx,nx),c(2:ny,ny)],
                          data[c(2:nx,nx),c(1,1:(ny-1))],
                          data[c(1,1:(nx-1)),c(2:ny,ny)],
                          data[c(1,1:(nx-1)),c(1,1:(ny-1))],
                          
                          data[c(3:nx,nx,nx),],
                          data[,c(3:ny,ny,ny)],
                          data[c(1,1,1:(nx-2)),],
                          data[,c(1,1,1:(ny-2))],
                          along = 3)
  
  dim(data_augmented) <- c(nx, ny, dim(data_augmented)[3], 1)
  
  
  ## predict missing values
  out <- Gapfill(data_augmented,
                 ## only predict missing values in first (original) image
                 subset = which(is.na(data_augmented[,,1,1])),
                 ## use parallel processing via R package foreach
                 dopar = TRUE,
                 ## tuning parameters of the algorithm
                 initialSize = c(2L, 2L, 100L, 100L),
                 nTargetImage = 2,
                 nQuant = 3,
                 ## restrict values to the following range
                 clipRange = range(data_augmented[,,1,1], na.rm=TRUE), 
                 ## return prediction interval
                 nPredict = 3, predictionInterval = TRUE  
  )
  
  ## extract prediction and prediction interval
  prediction <- out$fill[,,1,1,1]
  ciLo <- out$fill[,,1,1,2]; ciLo[is.na(ciLo)] <- prediction[is.na(ciLo)]
  ciUp <- out$fill[,,1,1,3]; ciUp[is.na(ciUp)] <- prediction[is.na(ciUp)]
  
  mean((Realdata >= matrix(ciLo,30,30)) * (Realdata <= matrix(ciUp,30,30)))
  
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  pred.Gapfill = prediction
  pred.NSBSR   = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  # pred.NSBSR   = matrix(as.vector(temp) + as.vector(temp2),30,30)
  
  summary((t(mylist[[k]])[1:30,1:30] - t(pred.Gapfill[1:30,1:30]) )^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.Gapfill[1:30,1:30])%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.Gapfill[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
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

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.Gapfill[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.Gapfill[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
int  = ciUp - ciLo; round(mean(int,na.rm=T),4)
time.Gapfill

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_Gapfill(toy).Rdata"))

