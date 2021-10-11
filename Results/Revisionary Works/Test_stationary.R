

# Generate a stationary image
# 
testimage <- matrix(rnorm(64*64), nrow=64, ncol=64)
fields::image.plot(testimage)
#
#Run test of stationarity

library(LS2Wstat)
TestofStat<-TOS2D(testimage)
TestofStat

library(geoR)

k <- 3
par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
brk = (seq(0.99*min(mylist[[k]],na.rm=T),1.01*max(mylist[[k]],na.rm=T),length.out=78))
x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[k]])),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(24,53,by=1),z=matrix(simdata$Y,30,30),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

sum(is.na(mylist[[k]]))/900 # [1] 0.1311111

mydata = data.frame(simdata$map,simdata$X,simdata$Y)
geod      = as.geodata(mydata,coords.col = 1:2, data.col = 10, covar.col = 7:9)
geod_full = as.geodata(mydata,coords.col = 1:2, data.col = 7, covar.col = 8:9)
fit <- likfit(geod,trend=trend.spatial(trend=~1,geodata=geod),cov.model="exponential",ini.cov.pars=c(1,1000),messages=FALSE)
estPara <- as.numeric(c(fit$beta,fit$nugget,fit$cov.pars))
krige = krige.control(type.krige = "ok",
                      trend.d = trend.spatial(trend=~1,geodata=geod),
                      trend.l = trend.spatial(trend=~1,geodata=geod_full),
                      cov.model = "exponential",nugget=fit$nugget,cov.pars=fit$cov.pars)
krige.UK = krige.conv(geodata=geod,locations=geod_full$coords,
                      krige=krige)
testY <- matrix(krige.UK$predict,30,30)

par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
brk = (seq(0.99*min(mylist[[k]],na.rm=T),1.01*max(mylist[[k]],na.rm=T),length.out=78))
x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
# fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[k]])),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
fields::image.plot(x=x,y=seq(24,53,by=1),z=matrix(testY,30,30),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

TestofStat1<-TOS2D(matrix(testY,30,30),nsamples=1000)
TestofStat1 # bootstrap p-value: 0.851

########################################################################################################################################

# load("C:/Users/user/Dropbox/200601StudyDataset.RData")

par(mfrow=c(1,1),mar=c(2,2,2,5))
colindex = rainbow(77,alpha=0.5)
x = seq(113,141,by=0.25);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n')
x = seq(112,142,by=0.25);y = 28/28*(x-113)+25
title(xlab="lon",ylab="lat",line=2)
fields::image.plot(x=x,y=y,z=matrix(O3L3.summer2019,(xrange+1),(yrange+1)),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
# fields::image.plot(x=x,y=seq(24,53,by=1),z=matrix(simdata$Y,30,30),col=colindex[6:68],axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

TestofStat2<-TOS2D(matrix(O3L3.summer2019,(xrange+1),(yrange+1)),nsamples=1000)
TestofStat2 # bootstrap p-value: 0.672
