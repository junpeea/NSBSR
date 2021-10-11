
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
fields::image.plot(O3L3.full,zlim=c(270,350))

#######################################################################################################################
# FRK
#######################################################################################################################

library(sp)
library(FRK)

time.FRK <- system.time({
  
  temp = expand.grid(c(x1:x2),c(y1:y2))
  mysimdata$x = temp$Var1;mysimdata$y = temp$Var2
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  coordinates(train.sim_data) = ~x + y # change into an sp object
  
  # ## Make BAUs as SpatialPixels
  
  sim_process <- mydata
  
  grid_BAUs <- sim_process                # assign BAUs
  
  grid_BAUs$Missing <- is.na(grid_BAUs$O3L3)      # mark which BAUs contain missing data
  
  grid_BAUs$O3L3  <- NULL                     # remove data from BAUs
  grid_BAUs$O3  <- NULL                     # remove data from BAUs
  grid_BAUs$NO2   <- NULL                     # remove data from BAUs
  grid_BAUs$VOC   <- NULL                     # remove data from BAUs
  grid_BAUs$R     <- NULL                     # remove data from BAUs
  # grid_BAUs$V     <- NULL                     # remove data from BAUs
  grid_BAUs$S     <- NULL                     # remove data from BAUs
  grid_BAUs$UV    <- NULL                     # remove data from BAUs
  # grid_BAUs$SOl   <- NULL                     # remove data from BAUs
  grid_BAUs$SO2   <- NULL                     # remove data from BAUs
  
  grid_BAUs$fs <- 1                          # set fs variation to unity
  
  coordinates(grid_BAUs)  <- ~x + y        # convert to SpatialPointsDataFrame
  
  gridded(grid_BAUs) <- TRUE                 # convert to SpatialPixelsDataFrame
  
  
  G <- auto_basis(manifold = plane(),
                  data=train.sim_data,
                  nres = 3,
                  regular = 1,
                  type = "bisquare",
                  subsamp = 20000)
  
  f <- z ~ 1
  
  S <- FRK(f = f,                       # formula for SRE model
           
           # data = sim_data,             # data
           data = train.sim_data,             # data
           
           basis = G,               # Basis
           
           BAUs = grid_BAUs,                 # BAUs
           
           tol = 1e-02)                   # EM iterations
  
  ### Check fit info
  
  ### Predict over BAUs
  grid_BAUs <- predict(S)
  
  ## Predict
  BAUs_pred <- SRE.predict(S)           # predict over all BAUs
  BAUs_pred_df <- data.frame(BAUs_pred) # convert to data frame
  
  ## Compute variance of predicted observations and conf. intervals
  BAUs_pred_df$sd_obs <- sqrt(BAUs_pred_df$sd^2 + S@Ve[1,1])
  int <- 2*1.96*BAUs_pred_df$sd_obs
  
  pred.FRK   = matrix(X.pred%*%fit0$coefficients + grid_BAUs$mu,(xrange+1),(yrange+1))
  
})

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(BAUs_pred_df$sd_obs,(xrange+1),(yrange+1))
mean((Realdata > pred.FRK - 1.96 * sdhatmat) * (Realdata < pred.FRK + 1.96 * sdhatmat))


par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.FRK,zlim=c(270,350))

mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.FRK))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.FRK))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.FRK) )^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.FRK) )^2,na.rm=T));round(rmse,3)
round(mean(int),3)
time.FRK

setwd("~/Desktop/Code3(Acting)")
# save.image(paste0(Sys.Date(),"_FRK.Rdata"))
# save.image(paste0(Sys.Date(),"_FRK2018.Rdata"))
save.image(paste0(Sys.Date(),"_FRK2017.Rdata"))

