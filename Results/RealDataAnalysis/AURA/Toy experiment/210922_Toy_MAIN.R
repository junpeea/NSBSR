########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Script name : NSBSR_Result4

# Date : 201908xx

# Author : YB JUN / CY Lim
########################################################################################################
########################################################################################################

############################################################################################## SETTINGS
rm(list = ls());
# setwd("C:/Users/user/Desktop/180102_LCY_THESIS") # set your main directory
setwd("D:/NSBSR_FINAL/Library")
source(file="PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="NSBSR_Func.R")             # Functions for Nonarametric Spectral Bayesian Spatial Regression
require(rgdal)
require(gdalUtils)
require(raster)
require(ggmap)
# image plot over googlemap
register_google(key="AIzaSyDg2IbZ7gGgQMBNtAxP5wBtgMLWiDHMIPg")
# register_google(key="AIzaSyCjtzlHwGhbI0qqwF5AukJ_LcjfZ0P-ffY")
# register_google(key="AIzaSyA0OzVQNRlyCDldlUSIHynNuomVT9REE0w")
mymap  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "satellite");plot(mymap)
mymap2  = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "terrain-labels");plot(mymap2)
mymap3 = get_map(location = c(lon = 127.024612, lat = 37.532600), zoom = 5, maptype = "toner-lines");plot(mymap3)


#######################################################################################################################
# Reference : NSBSR
# Reference : PBSR00
#######################################################################################################################

mse_r2_table = matrix(NA,9,4)
for(k in 1:9){
  # setwd("C:/Users/user/Desktop/NSBSR_FINAL/THESIS_Results(New)/Result4")
  setwd("D:/NSBSR_FINAL/THESIS_Results(New)/Result4")
  load(file=paste0("2019-08-06_MODISresult1_nd=",k,".Rdata"))
  Pred.summary = summary(mcmc.list(NSBSR.PredMat))
  n = nrow(Pred.summary$statistics)/3
  temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
  temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
  temp3 = temp + temp2
  
  sc = sqrt(mean(t(mylist[[k]])[1:30,1:30] - t(temp3[1:30,1:30]),na.rm=T)^2)
  r2 = summary(lm(as.vector(t(mylist[[k]])[1:30,1:30]) ~ as.vector(t(temp3))))$r.squared
  mse_r2_table[k,1]<-round(sc,4)
  mse_r2_table[k,3]<-round(r2,2)
  
  load(file=paste0("2019-08-06_MODISresult2_nd=",k,".Rdata"))
  Pred.summary = summary(mcmc.list(PBSR00.PredMat))
  
  n = nrow(Pred.summary$statistics)/3
  temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
  temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
  temp3 = temp + temp2
  
  sc = sqrt(mean(t(mylist[[k]])[1:30,1:30] - t(temp3[1:30,1:30]),na.rm=T)^2)
  r2 = summary(lm(as.vector(t(mylist[[k]])[1:30,1:30]) ~ as.vector(t(temp3))))$r.squared
  mse_r2_table[k,2]<-round(sc,4)
  mse_r2_table[k,4]<-round(r2,2)
}
mse_r2_table

require(coda)
setwd("D:/NSBSR_FINAL/THESIS_Results(New)/Result4")
load(file=paste0("2019-08-06_MODISresult1_nd=",1,".Rdata"))
Pred.summary = summary(mcmc.list(NSBSR.PredMat))
n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2
par(mfrow=c(2,2))
fields::image.plot(Y0[1:29,1:29])
fields::image.plot(Y00[1:29,1:29])
fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)])
fields::image.plot(temp3[1:29,1:29])


#######################################################################################################################
# Load data
#######################################################################################################################
load.image("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS/210922_MODIS_NSBSR_result.RData")
Y0  = matrix(Realdata,30,30)
Y00 = matrix(NA,30,30)
Y00[seq(1,30,by=2),seq(1,30,by=2)] = Y0[seq(1,30,by=2),seq(1,30,by=2)]
par(mfrow=c(2,2))
fields::image.plot(Y0[1:29,1:29])
fields::image.plot(Y00[1:29,1:29])
fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)])
fields::image.plot(Yhat[1:29,1:29])

