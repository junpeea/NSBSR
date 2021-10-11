
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
# Gapfill
#######################################################################################################################

## all packages are available on CRAN
library("gapfill") 
library("doParallel");
registerDoParallel(4) # run 4 tasks in parallel
library("abind")

time.Gapfill <- system.time({
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  sim_process = mydata
  sim_process$z = rep(NA,nrow(sim_process))
  sim_process$z[!is.na(mydata$O3L3)] <- train.sim_data$z
  # set.seed(1)
  # n.test = floor(nrow(sim_data)*0.2)
  # test.ind = sample(which(obs==TRUE),n.test)
  # sim_process2$z[test.ind] <- NA
  
  ## rearrange data as matrix
  data <- with(sim_process,
               array(z, c(length(unique(x)), length(unique(y)))))
  dim <- dim(data)
  nx <- dim[1]; ny <- dim[2]
  
  ## Imputation for Gapfill implementation
  temp = data[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)]
  temp[is.na(temp)] <- mean(data,na.rm=T)
  data[seq(1,xrange+1,by=2),seq(1,xrange+1,by=2)] <- temp
  
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
                 #subset = which(is.na(MissingValueHandling(O3L3.20190101))),
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
  ciLo <- out$fill[,,1,1,2]
  ciUp <- out$fill[,,1,1,3]
  
  pred.Gapfill = matrix(X.pred%*%fit0$coefficient+as.vector(prediction),nrow=(xrange+1),ncol=(yrange+1))
  
})

par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.Gapfill,zlim=c(270,350))

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = (ciUp-ciLo)/2; sdhatmat[is.na(sdhatmat)] = 0
mean((Realdata > pred.Gapfill - 1.96 * sdhatmat) * (Realdata < pred.Gapfill + 1.96 * sdhatmat))


mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.Gapfill))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.Gapfill))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.Gapfill) )^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.Gapfill) )^2,na.rm=T));round(rmse,3)
int  = ciUp - ciLo; round(mean(int,na.rm=T),3)
time.Gapfill

setwd("~/Desktop/Code3(Acting)")
save.image(paste0(Sys.Date(),"_Gapfill.Rdata"))

pdf(paste0(Sys.Date(),"_AURA_Gapfill.pdf"),width=7*1.5,height=7*1.5)
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
fields::image.plot(x=x,y=y,z=matrix(pred.Gapfill,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat",zlim=c(270,350))
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)
dev.off()

