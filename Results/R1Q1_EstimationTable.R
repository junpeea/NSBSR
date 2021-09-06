
rm(list = ls());
setwd("D:/NSBSR_FINAL")
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  
load(file="THESIS_Results(New)/Result1/190801_SPD_Iso32.Rdata")
# load(file="THESIS_Results(New)/Result1/190802_SPD_AnIso32.Rdata")
Esti.result = summary(mcmc.list(NSBSR.EstiMat))

ND = 100; D = 16; grid.size = c(32,32)
mspe    = matrix(NA,nrow=ND,ncol=7*6)
x = map.grid$x[which(simdata$map$obs==4)]
y = map.grid$y[which(simdata$map$obs==4)]
a = c(0,8,16)

thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)

resi    <- Y - as.matrix(X) %*% matrix(Esti.result$statistics[1:3,1],ncol=1)
sig2.e  <- as.numeric(var(resi))
sc.resi <- resi / sqrt(sig2.e)
e = matrix(sc.resi,ncol=33)
D = fftw2d(e); D2 = Re(D*Conj(D))
log.phi <- log(as.vector(D2/sum(D2)))
log.phi.mat = matrix(log.phi,33,33)
this.log.phi.mat = log.phi.mat[c(seq(1,32,by=2),1),c(seq(1,32,by=2),1)]

x = map.grid$x[which(simdata$map$obs==4)]
y = map.grid$y[which(simdata$map$obs==4)]
a = c(0,8,16)

pdf(paste0(Sys.Date(),"NSBSR-estimated-phi_IsoMatern05_Grid(32,32)_ND=",nd,"_3dlineci.pdf"))
# pdf(paste0(Sys.Date(),"NSBSR-estimated-phi_AnIsoMatern05_Grid(32,32)_ND=",nd,"_3dlineci.pdf"))

par(mfrow=c(1,1),mar=c(1,1,1,1))

# mylogdens2mat = logdens2mat[c(seq(1,32,by=2),1),c(seq(1,32,by=2),1)]
# lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(mylogdens2mat)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=F,type="o")
# lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(mylogdens2mat)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=T,type="o")
# lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(mylogdens2mat)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=T,type="o")

lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(this.log.phi.mat)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=1,colkey=FALSE,add=F,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(this.log.phi.mat)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=1,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(this.log.phi.mat)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=1,colkey=FALSE,add=T,type="o")

thetahat  = Esti.result$quantiles[(3+1+1):(3+1+ceiling(n^2/2)),1]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(mat.theta)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(mat.theta)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(mat.theta)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")

thetahat  = Esti.result$quantiles[(3+1+1):(3+1+ceiling(n^2/2)),3]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(mat.theta)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=3,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(mat.theta)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=3,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(mat.theta)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=3,colkey=FALSE,add=T,type="o")

thetahat  = Esti.result$quantiles[(3+1+1):(3+1+ceiling(n^2/2)),5]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(mat.theta)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(mat.theta)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(mat.theta)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")

dev.off()

# par(mfrow=c(1,1))
# lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(this.log.phi.mat)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=F,type="o")
# lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(this.log.phi.mat)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=T,type="o")
# lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(this.log.phi.mat)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=2,colkey=FALSE,add=T,type="o")
# lines3D(x=rep(a[1],length(which(x==a[1] & y<=16))),y=y[which(x==a[1] & y<=16)],z=as.vector(mat.theta)[which(x==a[1] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
# lines3D(x=rep(a[2],length(which(x==a[2] & y<=16))),y=y[which(x==a[2] & y<=16)],z=as.vector(mat.theta)[which(x==a[2] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")
# lines3D(x=rep(a[3],length(which(x==a[3] & y<=16))),y=y[which(x==a[3] & y<=16)],z=as.vector(mat.theta)[which(x==a[3] & y<=16)],xlim=c(0,16),ylim=c(0,16),zlim=c(-12,0),col=4,colkey=FALSE,add=T,type="o")


load(file="THESIS_Results(New)/Result1/190801_SPD_Iso64.Rdata")
# load(file="THESIS_Results(New)/Result1/190801_SPD_AnIso64.Rdata")
Esti.result = summary(mcmc.list(NSBSR.EstiMat))

thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)

resi    <- Y - as.matrix(X) %*% matrix(Esti.result$statistics[1:3,1],ncol=1)
sig2.e  <- as.numeric(var(resi))
sc.resi <- resi / sqrt(sig2.e)
e = matrix(sc.resi,ncol=65)
D = fftw2d(e); D2 = Re(D*Conj(D))
log.phi <- log(as.vector(D2/sum(D2)))
log.phi.mat = matrix(log.phi,65,65)
this.log.phi.mat = log.phi.mat[c(seq(1,64,by=2),1),c(seq(1,64,by=2),1)]

x = map.grid$x[which(simdata$map$obs==4)]
y = map.grid$y[which(simdata$map$obs==4)]
a = c(0,16,32)


pdf(paste0(Sys.Date(),"NSBSR-estimated-phi_IsoMatern05_Grid(64,64)_ND=",nd,"_3dlineci.pdf"))
# pdf(paste0(Sys.Date(),"NSBSR-estimated-phi_AnIsoMatern05_Grid(64,64)_ND=",nd,"_3dlineci.pdf"))

par(mfrow=c(1,1),mar=c(1,1,1,1))

lines3D(x=rep(a[1],length(which(x==a[1] & y<=32))),y=y[which(x==a[1] & y<=32)],z=as.vector(this.log.phi.mat)[which(x==a[1] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=1,colkey=FALSE,add=F,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=32))),y=y[which(x==a[2] & y<=32)],z=as.vector(this.log.phi.mat)[which(x==a[2] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=1,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=32))),y=y[which(x==a[3] & y<=32)],z=as.vector(this.log.phi.mat)[which(x==a[3] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=1,colkey=FALSE,add=T,type="o")

thetahat  = Esti.result$quantiles[(3+1+1):(3+1+ceiling(n^2/2)),1]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
lines3D(x=rep(a[1],length(which(x==a[1] & y<=32))),y=y[which(x==a[1] & y<=32)],z=as.vector(mat.theta)[which(x==a[1] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=32))),y=y[which(x==a[2] & y<=32)],z=as.vector(mat.theta)[which(x==a[2] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=32))),y=y[which(x==a[3] & y<=32)],z=as.vector(mat.theta)[which(x==a[3] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=4,colkey=FALSE,add=T,type="o")

thetahat  = Esti.result$quantiles[(3+1+1):(3+1+ceiling(n^2/2)),3]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
lines3D(x=rep(a[1],length(which(x==a[1] & y<=32))),y=y[which(x==a[1] & y<=32)],z=as.vector(mat.theta)[which(x==a[1] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=3,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=32))),y=y[which(x==a[2] & y<=32)],z=as.vector(mat.theta)[which(x==a[2] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=3,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=32))),y=y[which(x==a[3] & y<=32)],z=as.vector(mat.theta)[which(x==a[3] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=3,colkey=FALSE,add=T,type="o")

thetahat  = Esti.result$quantiles[(3+1+1):(3+1+ceiling(n^2/2)),5]
mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
lines3D(x=rep(a[1],length(which(x==a[1] & y<=32))),y=y[which(x==a[1] & y<=32)],z=as.vector(mat.theta)[which(x==a[1] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[2],length(which(x==a[2] & y<=32))),y=y[which(x==a[2] & y<=32)],z=as.vector(mat.theta)[which(x==a[2] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=4,colkey=FALSE,add=T,type="o")
lines3D(x=rep(a[3],length(which(x==a[3] & y<=32))),y=y[which(x==a[3] & y<=32)],z=as.vector(mat.theta)[which(x==a[3] & y<=32)],xlim=c(0,32),ylim=c(0,32),zlim=c(-15,0),col=4,colkey=FALSE,add=T,type="o")

dev.off()




beta.resultMat = matrix(NA,16*7,100)
for(nd in 1:100){
  
  setwd("D:/NSBSR_FINAL");ND <- 100
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result = summary(mcmc.list(NSBSR.EstiMat))
    beta.resultMat[c(1:3),nd] <- Esti.result$statistics[1:3,1]

  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Nugget_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
    Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
    Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
    Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
    beta.resultMat[c(4:6),nd]   <- Esti.result00$statistics[1:3,1]
    beta.resultMat[c(7:9),nd]   <- Esti.result01$statistics[1:3,1]
    beta.resultMat[c(10:12),nd] <- Esti.result05$statistics[1:3,1]
    beta.resultMat[c(13:15),nd] <- Esti.result20$statistics[1:3,1]
    
  }
  
  cat(nd,"1/7 Nugget procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result = summary(mcmc.list(NSBSR.EstiMat))
    beta.resultMat[16+c(1:3),nd] <- Esti.result$statistics[1:3,1]
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern01_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
    Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
    Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
    Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
    beta.resultMat[16+c(4:6),nd]   <- Esti.result00$statistics[1:3,1]
    beta.resultMat[16+c(7:9),nd]   <- Esti.result01$statistics[1:3,1]
    beta.resultMat[16+c(10:12),nd] <- Esti.result05$statistics[1:3,1]
    beta.resultMat[16+c(13:15),nd] <- Esti.result20$statistics[1:3,1]
  }
  
  cat(nd,"2/7 Matern01 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result = summary(mcmc.list(NSBSR.EstiMat))
    beta.resultMat[16*2+c(1:3),nd] <- Esti.result$statistics[1:3,1]
  }
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
    Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
    Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
    Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
    beta.resultMat[16*2+c(4:6),nd]   <- Esti.result00$statistics[1:3,1]
    beta.resultMat[16*2+c(7:9),nd]   <- Esti.result01$statistics[1:3,1]
    beta.resultMat[16*2+c(10:12),nd] <- Esti.result05$statistics[1:3,1]
    beta.resultMat[16*2+c(13:15),nd] <- Esti.result20$statistics[1:3,1]
  }
  
  cat(nd,"3/7 Matern05 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result = summary(mcmc.list(NSBSR.EstiMat))
    beta.resultMat[16*3+c(1:3),nd] <- Esti.result$statistics[1:3,1]
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp05_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
    Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
    Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
    Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
    beta.resultMat[16*3+c(4:6),nd]   <- Esti.result00$statistics[1:3,1]
    beta.resultMat[16*3+c(7:9),nd]   <- Esti.result01$statistics[1:3,1]
    beta.resultMat[16*3+c(10:12),nd] <- Esti.result05$statistics[1:3,1]
    beta.resultMat[16*3+c(13:15),nd] <- Esti.result20$statistics[1:3,1]
  }
  
  cat(nd,"4/7 SqExp05 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result = summary(mcmc.list(NSBSR.EstiMat))
    beta.resultMat[16*4+c(1:3),nd] <- Esti.result$statistics[1:3,1]
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_SqExp15_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
    Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
    Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
    Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
    beta.resultMat[16*4+c(4:6),nd]   <- Esti.result00$statistics[1:3,1]
    beta.resultMat[16*4+c(7:9),nd]   <- Esti.result01$statistics[1:3,1]
    beta.resultMat[16*4+c(10:12),nd] <- Esti.result05$statistics[1:3,1]
    beta.resultMat[16*4+c(13:15),nd] <- Esti.result20$statistics[1:3,1]
  }
  
  cat(nd,"5/7 SqExp15 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result = summary(mcmc.list(NSBSR.EstiMat))
    beta.resultMat[16*5+c(1:3),nd] <- Esti.result$statistics[1:3,1]
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Matern20_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
    Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
    Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
    Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
    beta.resultMat[16*5+c(4:6),nd]   <- Esti.result00$statistics[1:3,1]
    beta.resultMat[16*5+c(7:9),nd]   <- Esti.result01$statistics[1:3,1]
    beta.resultMat[16*5+c(10:12),nd] <- Esti.result05$statistics[1:3,1]
    beta.resultMat[16*5+c(13:15),nd] <- Esti.result20$statistics[1:3,1]
  }
  
  cat(nd,"6/7 Matern20 procedure complete!");print(Sys.time())
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/NSBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result = summary(mcmc.list(NSBSR.EstiMat))
    beta.resultMat[16*6+c(1:3),nd] <- Esti.result$statistics[1:3,1]
  }
  
  if(file.exists(paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))){
    load(file=paste0("THESIS_Results(New)/Result3/ISOresults/PBSR_Gauss_Grid(32,32)_ND=",ND,"_nd=",nd,".Rdata"))
    Esti.result00 = summary(mcmc.list(PBSR00.EstiMat))
    Esti.result01 = summary(mcmc.list(PBSR01.EstiMat))
    Esti.result05 = summary(mcmc.list(PBSR05.EstiMat))
    Esti.result20 = summary(mcmc.list(PBSR20.EstiMat))
    beta.resultMat[16*6+c(4:6),nd]   <- Esti.result00$statistics[1:3,1]
    beta.resultMat[16*6+c(7:9),nd]   <- Esti.result01$statistics[1:3,1]
    beta.resultMat[16*6+c(10:12),nd] <- Esti.result05$statistics[1:3,1]
    beta.resultMat[16*6+c(13:15),nd] <- Esti.result20$statistics[1:3,1]
  }
  
  cat(nd,"7/7 Gauss procedure complete!");print(Sys.time())
  
}
dim(beta.resultMat)
beta.resultMat
save(beta.resultMat,file=paste0(Sys.Date(),"_beta_resultMat.Rdata"))

beta.trueMat = matrix(rep(rep(c(rep(c(0.01,0.02,0.03),5),NA),7),100),ncol=100)
result.beta.mse.mat = (beta.resultMat - beta.trueMat)^2
result.beta.rmse.vec = sqrt(rowMeans(result.beta.mse.mat))
result.beta.Mat = matrix(NA,nrow=4*2,nco=6*4)
result.beta.Mat[1:3,1:5] <- matrix(result.beta.rmse.vec[c(1:15)],nrow=3)
result.beta.Mat[1:3,6+c(1:5)] <- matrix(result.beta.rmse.vec[16+c(1:15)],nrow=3)
result.beta.Mat[1:3,6*2+c(1:5)] <- matrix(result.beta.rmse.vec[16*2+c(1:15)],nrow=3)
result.beta.Mat[1:3,6*3+c(1:5)] <- matrix(result.beta.rmse.vec[16*3+c(1:15)],nrow=3)
result.beta.Mat[4+c(1:3),1:5] <- matrix(result.beta.rmse.vec[16*4+c(1:15)],nrow=3)
result.beta.Mat[4+c(1:3),6+c(1:5)] <- matrix(result.beta.rmse.vec[16*5+c(1:15)],nrow=3)
result.beta.Mat[4+c(1:3),6*2+c(1:5)] <- matrix(result.beta.rmse.vec[16*6+c(1:15)],nrow=3)
result.beta.Mat
save(result.beta.Mat,file=paste0(Sys.Date(),"_result_beta_Mat.Rdata"))
write.csv(round(result.beta.Mat,4),file=paste0(Sys.Date(),"_result_beta_rmse_Mat.csv"))


beta.resultVec = rowMeans(beta.resultMat,na.rm=T)
result.beta.Mat = matrix(NA,nrow=4*2,nco=6*4)
result.beta.Mat[1:3,1:5] <- matrix(beta.resultVec[c(1:15)],nrow=3)
result.beta.Mat[1:3,6+c(1:5)] <- matrix(beta.resultVec[16+c(1:15)],nrow=3)
result.beta.Mat[1:3,6*2+c(1:5)] <- matrix(beta.resultVec[16*2+c(1:15)],nrow=3)
result.beta.Mat[1:3,6*3+c(1:5)] <- matrix(beta.resultVec[16*3+c(1:15)],nrow=3)
result.beta.Mat[4+c(1:3),1:5] <- matrix(beta.resultVec[16*4+c(1:15)],nrow=3)
result.beta.Mat[4+c(1:3),6+c(1:5)] <- matrix(beta.resultVec[16*5+c(1:15)],nrow=3)
result.beta.Mat[4+c(1:3),6*2+c(1:5)] <- matrix(beta.resultVec[16*6+c(1:15)],nrow=3)
save(result.beta.Mat,file=paste0(Sys.Date(),"_result_beta_Mat.Rdata"))

true.beta.vec = rep(c(rep(c(0.01,0.02,0.03,NA),5),c(NA,NA,NA,NA)),8)
true.beta.mat = matrix(true.beta.vec,nrow=4)
true.beta.mat = rbind(true.beta.mat[,1:24],true.beta.mat[,25:48])
result.beta.bias.mat = round(result.beta.Mat-true.beta.mat,4)
write.csv(result.beta.bias.mat,file=paste0(Sys.Date(),"_result_beta_bias_Mat.csv"))



