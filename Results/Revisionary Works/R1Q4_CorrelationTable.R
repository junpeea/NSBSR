

# Experiment: Isotropic32
rm(list = ls());
setwd("D:/NSBSR_FINAL")
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  
library(spam);library(dplyr);library(fields)
load(file="THESIS_Results(New)/Result1/190801_SPD_Iso32.Rdata")
A = diag(2)
# load(file="THESIS_Results(New)/Result1/190802_SPD_AnIso32.Rdata")
# A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)


Esti.result = summary(mcmc.list(NSBSR.EstiMat))
thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
mat.lambda = exp(matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n))

corr = Re(fftw2d(mat.lambda,inverse=1)); corr[corr<0] <- 0
temp = toeplitz.spam(c(1:17))%>%as.matrix
temp2 = lapply(c(1:17),function(w) toeplitz.spam(corr[,w]))
cor_mtx = matrix(NA,289,289)
for(i in 1:289){
  for(j in 1:289){
    a = i%%17; b = ceiling(i/17)
    c = j%%17; d = ceiling(j/17)
    if(a==0) a=17; if(c==0) c=17;
    cor_mtx[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

gp2=gp(simdata$D,matern.specdens,c(simdata$true.cov.par[3]/simdata$D[1],simdata$true.cov.par[1]))
distmat = as.matrix(dist(expand.gpgrid(gp2)))

distvec = as.numeric(as.matrix(dist(as.matrix(gp2$omega)%*%A))[,1])
r = Matern.cor.to.range(16, nu=0.5, cor.target=0.001)
corr2mat = matrix(Matern(distvec,smoothness=0.5,range=r),ncol=32)


temp = toeplitz.spam(c(1:17))%>%as.matrix
temp2 = lapply(c(1:17),function(w) toeplitz.spam(corr2mat[,w]))
cor_mtx_true = matrix(NA,289,289)
for(i in 1:289){
  for(j in 1:289){
    a = i%%17; b = ceiling(i/17)
    c = j%%17; d = ceiling(j/17)
    if(a==0) a=17; if(c==0) c=17;
    cor_mtx_true[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

pdf(paste0(Sys.Date(),"_True_IsoMatern05_Grid(32,32)_ND=",nd,"_GammaMatrix.pdf"))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx_true,zlim=c(0,1))
dev.off()
pdf(paste0(Sys.Date(),"_NSBSR-estimated_IsoMatern05_Grid(32,32)_ND=",nd,"_GammaMatrix.pdf"))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx,zlim=c(0,1))
dev.off()


# Experiment: AnIsotropic32
rm(list = ls());
setwd("D:/NSBSR_FINAL")
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  
library(spam);library(dplyr);library(fields)
load(file="THESIS_Results(New)/Result1/190802_SPD_AnIso32.Rdata")
A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)


Esti.result = summary(mcmc.list(NSBSR.EstiMat))
thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
mat.lambda = exp(matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n))

corr = Re(fftw2d(mat.lambda,inverse=1)); corr[corr<0] <- 0
temp = toeplitz.spam(c(1:n))%>%as.matrix
temp2 = lapply(c(1:n),function(w) toeplitz.spam(corr[,w]))
cor_mtx = matrix(NA,n^2,n^2)
for(i in 1:n^2){
  for(j in 1:n^2){
    a = i%%n; b = ceiling(i/n)
    c = j%%n; d = ceiling(j/n)
    if(a==0) a=n; if(c==0) c=n;
    cor_mtx[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

gp2=gp(simdata$D,matern.specdens,c(simdata$true.cov.par[3]/simdata$D[1],simdata$true.cov.par[1]))
distmat = as.matrix(dist(expand.gpgrid(gp2)))

distvec = as.numeric(as.matrix(dist(as.matrix(gp2$omega)%*%A))[,1])
r = Matern.cor.to.range(16, nu=0.5, cor.target=0.05)
corr2mat = matrix(Matern(distvec,smoothness=0.5,range=r),ncol=32)


temp = toeplitz.spam(c(1:n))%>%as.matrix
temp2 = lapply(c(1:n),function(w) toeplitz.spam(corr2mat[,w]))
cor_mtx = matrix(NA,n^2,n^2)
for(i in 1:n^2){
  for(j in 1:n^2){
    a = i%%n; b = ceiling(i/n)
    c = j%%n; d = ceiling(j/n)
    if(a==0) a=n; if(c==0) c=n;
    cor_mtx[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

pdf(paste0(Sys.Date(),"_True_AnIsoMatern05_Grid(32,32)_ND=",nd,"_GammaMatrix.pdf"))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx_true,zlim=c(0,1))
dev.off()
pdf(paste0(Sys.Date(),"_NSBSR-estimated_AnIsoMatern05_Grid(32,32)_ND=",nd,"_GammaMatrix.pdf"))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx,zlim=c(0,1))
dev.off()


  
# Experiment: Isotropic64
rm(list = ls());
setwd("D:/NSBSR_FINAL")
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  
library(spam);library(dplyr);library(fields)
load(file="THESIS_Results(New)/Result1/190801_SPD_Iso64.Rdata")
A = diag(2)
# load(file="THESIS_Results(New)/Result1/190802_SPD_AnIso32.Rdata")
# A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)


Esti.result = summary(mcmc.list(NSBSR.EstiMat))
thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
mat.lambda = exp(matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n))

corr = Re(fftw2d(mat.lambda,inverse=1)); corr[corr<0] <- 0
temp = toeplitz.spam(c(1:n))%>%as.matrix
temp2 = lapply(c(1:n),function(w) toeplitz.spam(corr[,w]))
cor_mtx = matrix(NA,n^2,n^2)
for(i in 1:n^2){
  for(j in 1:n^2){
    a = i%%n; b = ceiling(i/n)
    c = j%%n; d = ceiling(j/n)
    if(a==0) a=n; if(c==0) c=n;
    cor_mtx[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

gp2=gp(simdata$D,matern.specdens,c(simdata$true.cov.par[3]/simdata$D[1],simdata$true.cov.par[1]))
distmat = as.matrix(dist(expand.gpgrid(gp2)))

distvec = as.numeric(as.matrix(dist(as.matrix(gp2$omega)%*%A))[,1])
r = Matern.cor.to.range(32, nu=0.5, cor.target=0.0001)
corr2mat = matrix(Matern(distvec,smoothness=0.5,range=r),ncol=64)

temp = toeplitz.spam(c(1:n))%>%as.matrix
temp2 = lapply(c(1:n),function(w) toeplitz.spam(corr2mat[,w]))
cor_mtx_true = matrix(NA,n^2,n^2)
for(i in 1:n^2){
  for(j in 1:n^2){
    a = i%%n; b = ceiling(i/n)
    c = j%%n; d = ceiling(j/n)
    if(a==0) a=n; if(c==0) c=n;
    cor_mtx_true[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

pdf(paste0(Sys.Date(),"_True_IsoMatern05_Grid(64,64)_ND=",nd,"_GammaMatrix.pdf"))
# par(mfrow=c(1,2),mar=c(1,1,1,1))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx_true,zlim=c(0,1))
dev.off()
pdf(paste0(Sys.Date(),"_NSBSR-estimated_IsoMatern05_Grid(64,64)_ND=",nd,"_GammaMatrix.pdf"))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx,zlim=c(0,1))
dev.off()



# Experiment: AnIsotropic64
rm(list = ls());
setwd("D:/NSBSR_FINAL")
source(file="THESIS_Results(New)/Result3/PBSR_Func_MCMC.R")         # Functions for Parametric Bayesian Spatial Regression
source(file="THESIS_Results(New)/Result3/NSBSR_Func.R")  
library(spam);library(dplyr);library(fields)
load(file="THESIS_Results(New)/Result1/190801_SPD_AnIso64.Rdata")
A = matrix(c(-sqrt(5),0,+sqrt(5),1),2,2)


Esti.result = summary(mcmc.list(NSBSR.EstiMat))
thetahat  = Esti.result$statistics[(3+1+1):(3+1+ceiling(n^2/2))]
mat.lambda = exp(matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n))

corr = Re(fftw2d(mat.lambda,inverse=1)); corr[corr<0] <- 0
temp = toeplitz.spam(c(1:n))%>%as.matrix
temp2 = lapply(c(1:n),function(w) toeplitz.spam(corr[,w]))
cor_mtx = matrix(NA,n^2,n^2)
for(i in 1:n^2){
  for(j in 1:n^2){
    a = i%%n; b = ceiling(i/n)
    c = j%%n; d = ceiling(j/n)
    if(a==0) a=n; if(c==0) c=n;
    cor_mtx[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

gp2=gp(simdata$D,matern.specdens,c(simdata$true.cov.par[3]/simdata$D[1],simdata$true.cov.par[1]))
distmat = as.matrix(dist(expand.gpgrid(gp2)))

distvec = as.numeric(as.matrix(dist(as.matrix(gp2$omega)%*%A))[,1])
r = Matern.cor.to.range(32, nu=0.5, cor.target=0.075)
corr2mat = matrix(Matern(distvec,smoothness=0.5,range=r),ncol=64)


temp = toeplitz.spam(c(1:n))%>%as.matrix
temp2 = lapply(c(1:n),function(w) toeplitz.spam(corr2mat[,w]))
cor_mtx_true = matrix(NA,n^2,n^2)
for(i in 1:n^2){
  for(j in 1:n^2){
    a = i%%n; b = ceiling(i/n)
    c = j%%n; d = ceiling(j/n)
    if(a==0) a=n; if(c==0) c=n;
    cor_mtx_true[i,j] = temp2[[temp[a,c]]][b,d]
  }
}

pdf(paste0(Sys.Date(),"_True_AnIsoMatern05_Grid(32,32)_ND=",nd,"_GammaMatrix.pdf"))
# par(mfrow=c(1,2),mar=c(1,1,1,1))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx_true,zlim=c(0,1))
dev.off()
pdf(paste0(Sys.Date(),"_NSBSR-estimated_AnIsoMatern05_Grid(32,32)_ND=",nd,"_GammaMatrix.pdf"))
par(mfrow=c(1,1),mar=c(1,1,1,1))
fields::image.plot(cor_mtx,zlim=c(0,1))
dev.off()
