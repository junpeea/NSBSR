########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : Show_MCMC_NSBSR_Images.R

########################################################################################################
########################################################################################################
require(lattice)
library("plot3D")
library(spectralGP)
library(dplyr)


n1=simdata$D[1]; n2=simdata$D[2]; alpha=2.0
# gp2=gp(c(n1,n2),matern.specdens,c(0.5,alpha))
# dens2=matern.specdens(gp2$omega,c(0.5,alpha),d=2)
# gp2=gp(c(n1,n2),matern.specdens,c(1.2,alpha))
# dens2=matern.specdens(gp2$omega,c(1.2,alpha),d=2)
if(mode=="AnISO"){
   A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)  
  # A = matrix(c(1,0,0,0.5),2,2)  
}else{
  A = diag(2)
}
gp2=gp(c(n1,n2),matern.specdens,c(0.2,alpha))
dens2=matern.specdens(as.matrix(gp2$omega)%*%A,c(0.25,alpha),d=2)
dens2mat = matrix(dens2,n1,n2)
logdens2mat = log(dens2mat / sum(dens2mat))
nh = ceiling(n1/2)
re.logdens2mat = rbind(cbind(logdens2mat[rev(1:nh),rev(1:nh)],logdens2mat[rev(1:nh),1:nh]),
                       cbind(logdens2mat[1:nh,rev(1:nh)],logdens2mat[1:nh,1:nh]))
logdens2mat = rbind(cbind(re.logdens2mat[rev(1:nh),rev(1:nh)],re.logdens2mat[rev(1:nh),1:nh]),
                    cbind(re.logdens2mat[1:nh,rev(1:nh)],re.logdens2mat[1:nh,1:nh]))
#######################################################################
par(mfrow=c(2,3),mar=c(2,2,3,3))
X = simdata$X; Y = simdata$Y

Yhat = PredResult$result1 + PredResult$result2
brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
score = sqrt(mean((Y-Yhat)^2))

Ymat  = matrix(Y,   ncol=sqrt(simdata$nloc))
Yhmat = matrix(Yhat,ncol=sqrt(simdata$nloc))

col <- rev(rainbow(77)[1:63])

# Take 1
fields::image.plot(Ymat,main="TRUE",axes=F)
# Take 2
fields::image.plot(Yhmat,main=paste0("MSE=",round(score,2)),axes=F)
# Take 3
plot(Y,Yhat);abline(b=1,lty=2,col=2,axes=F)
# Take 4
# fields::image.plot(t(logdens2mat[c(seq(1,16,by=2),seq(18,32,by=2)),c(seq(1,16,by=2),seq(18,32,by=2))]),axes=F,cex.main=3.2,legend.mar=5.0,zlim=c(-13,-0.5))
# fields::image.plot(t(logdens2mat[seq(2,32,by=2),seq(2,32,by=2)]),axes=F,cex.main=3.2,legend.mar=5.0,zlim=c(-13,-0.5))
fields::image.plot(t(re.logdens2mat),axes=F,cex.main=3.2,legend.mar=5.0,zlim=c(-16,-0.5))
# fields::image.plot(re.logdens2mat[9:24,9:24],axes=F,cex.main=3.2,legend.mar=5.0,zlim=c(-11,-0.5))
# fields::image.plot(re.logdens2mat[seq(1,32,by=2),seq(1,32,by=2)],axes=F,cex.main=3.2,legend.mar=5.0,zlim=c(-11,-0.5))
# fields::image.plot(exp(re.logdens2mat),axes=F,cex.main=3.2,legend.mar=5.0,zlim=c(0,0.16))

# Take 5
meantheta = Esti.summary$statistics[(new.mcmcobj$n.beta+1+1):(new.mcmcobj$n.beta+1+new.mcmcobj$n.theta),1]
n = sqrt(simdata$nobs)
if(n%%2==0){
  mat.theta  = matrix(c(meantheta,rev(meantheta)),ncol=n)
}else{
  mat.theta  = matrix(c(meantheta,rev(meantheta[1:(length(meantheta)-1)])),ncol=n)
} 
nh = ceiling(n/2)

rearranged.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
                         cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))

# rearranged.theta = rbind(cbind(mat.theta[rev(c(8:16)),rev(c(8:16))],mat.theta[rev(c(8:16)),1:nh]),
#                          cbind(mat.theta[1:nh,rev(c(8:16))],mat.theta[1:nh,1:nh]))


temp = rearranged.theta
temp[1:3,7:10] <- rearranged.theta[2,2]
temp[14:16,7:10] <- rearranged.theta[2,2]
temp[7:10,1:3] <- rearranged.theta[2,2]
temp[7:10,14:16] <-rearranged.theta[2,2]
# fields::image.plot(temp,main="SPD(averaged)",legend.mar=5.0,axes=F,zlim=c(-14,-0.5))
# fields::image.plot(mat.theta,main="SPD(averaged)",legend.mar=5.0,axes=F,zlim=c(-12,-0.5))
fields::image.plot((rearranged.theta),main="SPD(averaged)",legend.mar=5.0,axes=F,zlim=c(-16,-0.5))
# fields::image.plot((rearranged.theta),main="SPD(averaged)",legend.mar=5.0,axes=F)
# fields::image.plot(exp(rearranged.theta),main="SPD(averaged)",legend.mar=5.0,axes=F,zlim=c(0,0.16))

# summary(1.35*(rearranged.theta+1.35)%>%as.vector)
# summary(re.logdens2mat[seq(1,32,by=2),seq(1,32,by=2)]%>%as.vector)

# Take 6
lambda = new.mcmcobj$Parset$lambda
mat.lambda = matrix(lambda,ncol=n)/sum(lambda)
# fields::image.plot(mat.lambda,main="exp(SPD)",axes=F)
# plot((rearranged.theta)%>%as.vector,type='o',ylim=c(-13,-1.5))
# plot((rearranged.theta*1.8+3.9)%>%as.vector,type='o',ylim=c(-11.5,-1.5))
# plot((rearranged.theta*90+492)%>%as.vector,type='o',ylim=c(-11.5,-1.5))
# plot(rearranged.theta%>%as.vector,type='o',ylim=c(-11.5,-1.5))
# par(new=T);plot(re.logdens2mat[seq(2,32,by=2),seq(2,32,by=2)]%>%as.vector,type='o',col=2,ylim=c(-13,-1.5))


wireframe(t(re.logdens2mat[c(seq(1,16,by=2),seq(18,32,by=2)),c(seq(1,16,by=2),seq(18,32,by=2))]),shade=TRUE,zlab="")
wireframe(rearranged.theta,shade=TRUE,zlab="")
wireframe(temp,shade=TRUE,zlab="")
# summary((rearranged.theta*1.8+3.9)%>%as.vector)
# summary(re.logdens2mat[seq(1,32,by=2),seq(1,32,by=2)]%>%as.vector)

x = map.grid$x[which(simdata$map$obs==4)]
y = map.grid$y[which(simdata$map$obs==4)]
a = c(0,8,16)
mylogdens2mat = t(logdens2mat[c(seq(1,16,by=2),seq(18,32,by=2)),c(seq(1,16,by=2),seq(18,32,by=2))])
lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(mylogdens2mat)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=F,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(mylogdens2mat)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(mylogdens2mat)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(mat.theta)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(mat.theta)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(mat.theta)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")



#######################################################################
rm(bdd); rm(obs); rm(D)


