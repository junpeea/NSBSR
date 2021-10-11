#######################################################################################################################
# Periodic Embedding
#######################################################################################################################

# devtools::install_github("joeguinness/aldodevel")
# devtools::install_github("joeguinness/npspec")
library(npspec)
setwd("~/GitHub/npspec/R")
source("iterative_estimation.R")
source("parametricfuns.R")
source("RcppExports.R")
library(dplyr)

setwd("~/GitHub/npspec/Example")
# open: "200206_Toy_Periodic.R"
load(file="200206_Toydata.Rdata")

time.periodic <- system.time({
  
  n <- 30
  temp = expand.grid(c(0:(n-1)),c(0:(n-1)))
  sim_process = data.frame(x=temp$Var1,y=temp$Var2)
  sim_process$Y = as.vector(Y0)
  ind = is.na(as.vector(Y00));summary(ind)
  sim_data <- sim_process[which(ind==FALSE),]
  sim_data$x1 = rep(1,nrow(sim_data))
  sim_data$x2 = scale(sim_process$x[which(ind==FALSE)])
  sim_data$x3 = scale(sim_process$y[which(ind==FALSE)])
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  sim_data$z = fit0$residuals
  sim_data$std <- 0.1
  
  Ymat   = matrix(Y00,30,30)
  indmat = !is.na(Ymat)
  Xmat   = matrix(NA,900,3)
  Xmat[,1] <-rep(1,nrow(sim_process))
  Xmat[,2] <-scale(sim_process$x)
  Xmat[,3] <-scale(sim_process$y)
  # Xmat[indmat,1] <-sim_data$x1
  # Xmat[indmat,2] <-sim_data$x2
  # Xmat[indmat,3] <-sim_data$x3
  Xarray = array(NA,dim=c(30,30,3))
  Xarray[,,1] <- matrix(Xmat[,1],30,30)
  Xarray[,,2] <- matrix(Xmat[,2],30,30)
  Xarray[,,3] <- matrix(Xmat[,3],30,30)
  args(iterate_spec)
  # fit.PeriodicEmbedd = iterate_spec(y=Ymat,observed=indmat,embed_fac = 2)
  fit.PeriodicEmbedd = iterate_spec(y=Ymat,X=Xarray,observed=indmat)
  pred.Periodic = matrix(fit.PeriodicEmbedd$condsim,30,30)
  
  cond_diff <- array(NA, dim(fit.PeriodicEmbedd$condsim) )
  for(j in 1:dim(fit.PeriodicEmbedd$condsim)[3]) cond_diff[,,j] <- fit.PeriodicEmbedd$condsim[,,j] - fit.PeriodicEmbedd$condexp
  meansq <- function(x) 1/length(x)*sum(x^2)
  predvar_mat <- apply( cond_diff, c(1,2), meansq )
  
  par(mfrow=c(1,2))
  fields::image.plot(Realdata)
  fields::image.plot(pred.Periodic)
  summary( predvar_mat)
  
  sdhatmat = sqrt(predvar_mat)
  mean((Realdata >= pred.Periodic - 1.96 * sdhatmat) * (Realdata <= pred.Periodic + 1.96 * sdhatmat))
  
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  
  # pred.NSBSR    = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  
  
})

summary((t(mylist[[k]])[1:30,1:30] - t(pred.NSBSR[1:30,1:30]))^2 %>% as.vector)
cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.NSBSR[1:30,1:30])%>%as.vector,use="complete.obs")
par(mfrow=c(2,2))
fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
fields::image.plot(pred.NSBSR[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.Periodic[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.Periodic[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int),4)
time.periodic

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_Periodic(toy).Rdata"))
