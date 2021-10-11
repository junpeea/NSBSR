########################################################################################################
########################################################################################################
# Thesis : Bayesian Spatial Process Regression with Nonparametric modeling of spectral densities
# Project name : NSBSR_Result4

# Date : 201811xx
# Author : YB JUN / Chae Young Lim
########################################################################################################
########################################################################################################

rm(list = ls());
load(file="THESIS_Results/Result4/CrossValid.Rdata")
pdf(paste0("THESIS_Results/Result4/",Sys.Date(),"MODIS_FIG1.pdf"),width=7*1.5,height=7*1.5)

colindex = rainbow(77,alpha=0.5)
brk = (seq(-4.5,4,length.out=78))
par(mfrow=c(2,2))

x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n',xlim=c(112.5,140.5),ylim=c(25,52))
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[1]])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

x = seq(112,141,by=2);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="TRAINING DATASET",xaxt='n',yaxt='n')
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=2),z=(t(mylist[[1]][seq(1,30,by=2),seq(1,30,by=2)])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


Pred.summary = summary(mcmc.list(NSBSR.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2


x = seq(112,141,by=1);y = 28/28*(x-113)+25
sc = mean(t(mylist[[1]])[2:30,2:30] - t(temp3),na.rm=T)^2
r2 = summary(lm(as.vector(t(mylist[[1]])[2:30,2:30]) ~ as.vector(t(temp3))))$r.squared
plot(x,y,type="n",xlab="",ylab="",main=paste0("ESTIMATED (NSBSR) OZONE"),xaxt='n',yaxt='n',xlim=c(113,140),ylim=c(25.5,52.5))
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(25,53,by=1),z=(t(temp3)),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

Pred.summary = summary(mcmc.list(PBSR.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$statistics[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$statistics[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2

x = seq(112,141,by=1);y = 28/28*(x-113)+25
sc = mean(t(mylist[[1]])[2:30,2:30] - t(temp3),na.rm=T)^2
r2 = summary(lm(as.vector(t(mylist[[1]])[2:30,2:30]) ~ as.vector(t(temp3))))$r.squared
plot(x,y,type="n",xlab="",ylab="",main=paste0("ESTIMATED (PBSR) OZONE"),xaxt='n',yaxt='n',xlim=c(113,140),ylim=c(25.5,52.5))
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(25,53,by=1),z=(t(temp3)),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

dev.off()






rm(list = ls());
load(file="THESIS_Results/Result4/Pred2Higher.Rdata")
pdf(paste0("THESIS_Results/Result4/",Sys.Date(),"MODIS_FIG2.pdf"),width=7*1.5,height=7*1.5)

par(mfrow=c(2,2))
colindex = rainbow(77,alpha=0.5)
brk = (seq(-4.5,4,length.out=78))

x = seq(112,141,by=1);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="ORIGINAL DATASET",xaxt='n',yaxt='n',xlim=c(113,140),ylim=c(25.5,52.5))
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=1),z=(t(mylist[[1]])),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


Pred.summary = summary(mcmc.list(NSBSR.PredMat))

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,3],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),3],ncol=sqrt(n))
temp3 = temp + temp2



x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE",xaxt='n',yaxt='n',xlim=c(113,140),ylim=c(25.5,52.5))
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")

par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,1],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),1],ncol=sqrt(n))
temp3 = temp + temp2


x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE (2.5% quantile)",xaxt='n',yaxt='n',xlim=c(113,140),ylim=c(25.5,52.5))
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)

n = nrow(Pred.summary$statistics)/3
temp = matrix(Pred.summary$quantiles[1:n,5],ncol=sqrt(n))
temp2 = matrix(Pred.summary$quantiles[(n+1):(2*n),5],ncol=sqrt(n))
temp3 = temp + temp2

x = seq(112,141,by=0.5);y = 28/28*(x-113)+25
plot(x,y,type="n",xlab="",ylab="",main="PREDICTED (NSBSR) OZONE (97.5% quantile)",xaxt='n',yaxt='n',xlim=c(113,140),ylim=c(25.5,52.5))
title(xlab="lon",ylab="lat",line=2)
axis(side=1,at=seq(115,141,by=5),labels=seq(115,141,by=5),line=-0.5,tick=FALSE)
axis(side=2,at=seq(25,53,by=5),labels=seq(25,53,by=5),line=-0.5,tick=FALSE)
fields::image.plot(x=x,y=seq(24,53,by=0.5),z=t(temp3),col=colindex,breaks=brk,axes = FALSE,add=T,xlab="lon",ylab="lat")
par(new=T);plot(mymap2)
par(new=T);plot(mymap3)


dev.off()


