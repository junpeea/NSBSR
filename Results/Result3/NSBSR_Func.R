########################################################################################################
########################################################################################################
# Thesis : Bayesian Spatial Process Regression with Nonparametric modeling of autocorrelation
# Project name : NSBSR_func

# Date : 201805xx
# Author : YB JUN / Chae Young Lim
########################################################################################################
########################################################################################################

require(spectralGP)
require(RandomFields)
require(geoR)
require(coda)
require(fftwtools)

########################################################################################################
### CONTENTS

### fftw2d(e,inverse=0) : FFT-solver
### Get2dgridmap(seed,lgs=32,resol=2) : Create a new field object
### SetMCMC(seed,nd,n.chain=3,n.iter,n.burn) : Create a mcmc object
### BlockToBCCB(blockmtxlist,indmtx) : BLOCK list -> BCCB matrix
### Extend.lambda(mat.lambda) : Interior lambda matrix -> Bounded lambda matrix
### Expand.iso(mat) : Expansion function with respect to an isotrpoy assumption

### Simulate(map,seed,nd,beta=c(1,2,3),sig2_eps=1,cor.par,model.name,model.option)
########################################### Create a new simulation data object
# Simulate.Spherical = function(map,seed,nd,beta=c(1,2,3),sig2_eps=1,cor.par)
# Simulate.SqExp     = function(map,seed,nd,beta=c(1,2,3),sig2_eps=1,cor.par)
# Simulate.Matern    = function(map,seed,nd,beta=c(1,2,3),sig2_eps=1,cor.par)
# Simulate.Gaussian  = function(map,seed,nd,beta=c(1,2,3),sig2_eps=1,cor.par)
###############################################################################  

### EvalTrueMSPE(ND,simdata) : Evaluate True MSPE
### CalcMSPE.BSR(map,data0,BSR.PredMat) : Calcluate (Emperical/Model) MSPE for BSR-method
### CalcMSPE.UK(map,geodata0,geodata00,alpha) : Calcluate (Emperical/Model) MSPE for UK-method

############################################################################################## CONTENTS
########################################################################################################

fftw2d = function(e,inverse=0){
  return(t(mvfft(t(mvfft(e,inverse)))))
}

Get2dgridmap = function(seed=1,lgs=32,resol=2){
  out = list()
  
  delta = 1/resol
  sp1 = seq(0,lgs,by=delta);sp2 = seq(0,lgs,by=delta)
  omega1 = (-pi+(2*pi)/lgs*sp1);omega2 = (-pi+(2*pi)/lgs*sp2)
  space.grid = expand.grid(sp1,sp2)
  omega.grid = expand.grid(omega1,omega2)
  
  nloc = nrow(space.grid)
  X0 = rep(1,nloc)
  set.seed(seed); p = rbinom(nloc,1,0.5)
  seed=seed+123
  set.seed(seed); X1 = p*rnorm(nloc,0,1) + (1-p)*rnorm(nloc,0,5)
  seed=seed+123
  set.seed(seed); X2 = rexp(nloc,1)
  
  # bddind : boundary coordinates
  bddind = rep(0,nloc)
  bddind[which(omega.grid[,1]==0 | abs(omega.grid[,1])==pi | omega.grid[,2]==0 | abs(omega.grid[,2])==pi)] <- 1
  # bddind2 : auxilary boundary corrdinates
  bddind2 = rep(0,nloc)
  bddind2[which(space.grid[,1]==(lgs-0) | space.grid[,2]==(lgs-0))] <- 1
  bddind2[which(space.grid[,1]==(lgs-delta) | space.grid[,2]==(lgs-delta))] <- 1
  
  # obsind : observable coordinates
  obsind = rep(NA,nloc)
  obsind[ which( (space.grid[,1]*resol)%%resol==0 & (space.grid[,2]*resol)%%resol==0 ) ] <- 1
  obsind[ which(bddind2==1) ] <- NA
  # obsind2 : observed coordinates
  obsind2 = obsind
  seed=seed+123
  set.seed(seed); misind = rbinom(nloc,1,prob=0); obsind2[which(obsind==1 & misind==1)] <- 0
  
  fieldgrid = data.frame(cbind(space.grid,omega.grid,obsind,obsind2,bddind,bddind2,X0,X1,X2))
  colnames(fieldgrid) = c("x","y","w.x","w.y","obs","obs2","bdd","bdd2","X0","X1","X2")
  
  # par(mfrow=c(2,2))
  # plot(x=fieldgrid$x,y=fieldgrid$y,col=(fieldgrid$bdd+1),pch=20)
  # plot(x=fieldgrid$x,y=fieldgrid$y,col=(fieldgrid$bdd2+1),pch=20)
  # plot(x=fieldgrid$x,y=fieldgrid$y,col=(fieldgrid$obs+2),pch=20)
  # plot(x=fieldgrid$x,y=fieldgrid$y,col=(fieldgrid$obs2+2),pch=20)
  
  dist.obs.omega = as.matrix(dist(fieldgrid[which(fieldgrid$bdd2==0 & fieldgrid$obs==1),c("w.x","w.y")]))
  dist.obs.omega.x = as.matrix(dist(fieldgrid[which(fieldgrid$bdd2==0 & fieldgrid$obs==1),c("w.x")]))
  dist.obs.omega.y = as.matrix(dist(fieldgrid[which(fieldgrid$bdd2==0 & fieldgrid$obs==1),c("w.y")]))
  
  out$delta = delta
  out$gridsize = lgs
  out$nloc = nloc
  out$nobs = length(which(obsind==1))
  out$npred= length(which(bddind2==0))
  out$fieldgrid  = fieldgrid
  out$d.o.o = dist.obs.omega
  out$d.o.o.x = dist.obs.omega.x
  out$d.o.o.y = dist.obs.omega.y
  return(out)
}

SetMCMC = function(seed,nd,n.chain=3,n.iter,n.burn){
  out = list()
  for(chain in 1:n.chain){
    out[[chain]] = list()
    out[[chain]]$n.iter  = n.iter
    out[[chain]]$n.burn  = n.burn
    out[[chain]]$seed    = (nd*10+chain)*seed
    out[[chain]]$EstiMat = mcmc.list()
    out[[chain]]$PredMat = mcmc.list()
  }
  return(out)
}

BlockToBCCB = function(blockmtxlist,indmtx){
  n = ncol(indmtx)
  first.block <- blockmtxlist[[indmtx[1,1]]]
  for(col in 2:n){
    first.block = cbind(first.block,blockmtxlist[[indmtx[1,col]]])
  }
  finalmtx <- first.block
  for(row in 2:n){
    current.block <- blockmtxlist[[indmtx[row,1]]]
    for(col in 2:n){
      current.block = cbind(current.block,blockmtxlist[[indmtx[row,col]]])
    }
    finalmtx = rbind(finalmtx,current.block)
  }
  return(finalmtx)
}

Extend.lambda = function(mat.lambda){
  clambda = cbind(mat.lambda,mat.lambda[,1])
  lambda0 = rbind(clambda,clambda[1,])
  return(lambda0)
}

Expand.iso = function(mat){
  a = ncol(mat)-1
  rmat = rbind(mat[1:a,],apply(mat,2,function(w) rev(w)))
  cmat = cbind(rmat[,1:a],t(apply(t(rmat),2,function(w) rev(w))))
  return(cmat)
}

Simulate = function(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,model.name,model.option){
  
  out = list()
  
  delta = map$delta
  gridsize = map$gridsize
  nloc = map$nloc
  nobs = map$nobs
  fieldgrid = map$fieldgrid
  X = as.matrix(cbind(map$fieldgrid$X0,map$fieldgrid$X1,map$fieldgrid$X2))
  
  epsi = rep(NA,nloc)
  set.seed(seed*nd)
  epsi = RandomFields::RFsimulate(
    model=model.option,
    x=fieldgrid$x,y=fieldgrid$y
  )$variable1
  
  Y = X%*%beta + sqrt(sig2_eps) * epsi
  
  out$map <- map
  out$beta = beta
  out$cor.model = model.name
  out$true.cov.pars = c(cor.par[1],sig2_eps,cor.par[2])
  out$n.beta  = length(beta)
  out$n.epsi  = map$npred
  out$nd = nd
  out$DB = data.frame(cbind(fieldgrid,epsi,Y))
  return(out)
  
}

Simulate.Spheric = function(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,aniso.option=RMangle(angle=0,ratio=1)){
  
  model.name   = "Spheric"
  model.option = RMspheric(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  return(Simulate(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,model.name,model.option))
  
}

Simulate.SqExp   = function(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,aniso.option=RMangle(angle=0,ratio=1)){
  
  model.name   = "SquaredExponential"
  model.option = RMpoweredexp(alpha=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  return(Simulate(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,model.name,model.option))
  
}

Simulate.Matern    = function(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,aniso.option=RMangle(angle=0,ratio=1)){
  
  model.name   = "Matern"
  model.option = RMmatern(nu=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  return(Simulate(map,seed,nd,beta=beta,sig2_eps=sig2_eps,cor.par,model.name,model.option))
  
}

Simulate.Gauss = function(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,aniso.option=RMangle(angle=0,ratio=1)){
  
  model.name   = "Gauss"
  model.option = RMgauss(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  return(Simulate(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,model.name,model.option))
  
}

Simulate.Nugget = function(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,aniso.option=RMangle(angle=0,ratio=1)){
  
  model.name   = "Nugget"
  model.option = RMnugget(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  return(Simulate(map,seed,nd,beta=c(0.01,0.02,0.03),sig2_eps=1,cor.par,model.name,model.option))
  
}

EvalTrueMSPE = function(ND,simdata){
  
  out = rep(0,8)
  
  data0     = simdata$DB;                            
  data00    = simdata$DB[which(simdata$DB$obs==1),]; 
  dist1 = as.matrix(dist(data0[,1:2]))
  dist2 = as.matrix(dist(data00[,1:2]))
  
  cormtx1 = cov.spatial(dist1,cov.model="pure.nugget",cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="pure.nugget",cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  out[1] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  
  cormtx1 = cov.spatial(dist1,cov.model="spherical",cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="spherical",cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  out[2] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  
  cormtx1 = cov.spatial(dist1,cov.model="matern",kappa=0.1,cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="matern",kappa=0.1,cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  out[3] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  
  cormtx1 = cov.spatial(dist1,cov.model="matern",kappa=0.5,cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="matern",kappa=0.5,cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  out[4] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  
  cormtx1 = cov.spatial(dist1,cov.model="matern",kappa=2.0,cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="matern",kappa=2.0,cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  out[5] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  
  cormtx1 = cov.spatial(dist1,cov.model="powered.exponential",kappa=0.5,cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="powered.exponential",kappa=0.5,cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  out[6] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  
  cormtx1 = cov.spatial(dist1,cov.model="powered.exponential",kappa=1.5,cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="powered.exponential",kappa=1.5,cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  out[7] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  
  cormtx1 = cov.spatial(dist1,cov.model="gaussian",cov.pars=c(1,10))
  cormtx2 = cov.spatial(dist2,cov.model="gaussian",cov.pars=c(1,10))
  cormtx3 = cormtx1[,which(simdata$DB$obs==1)]
  if(is.singular.matrix(cormtx2)==FALSE){
    out[8] = mean( diag( cormtx1 - cormtx3 %*% solve(cormtx2)%*%t(cormtx3) ) )
  }
  return(out)
}

CalcMSPE.BSR = function(map,data0,BSR.PredMat){
  
  out = c(NA,NA)
  
  Y = matrix(data0$Y,ncol=sqrt(map$nloc))
  Ytrue = Y[3:nrow(Y),3:ncol(Y)]
  
  Pred.summary = summary(mcmc.list(BSR.PredMat))
  
  muhat = Pred.summary$statistics[1:(map$npred),1]
  ehat  = Pred.summary$statistics[(map$npred+1):(map$npred+map$npred),1]
  evar  = Pred.summary$statistics[(map$npred+map$npred+1):(map$npred+map$npred+map$npred),1]
  
  muhat = matrix(muhat,ncol=sqrt(length(muhat)))
  ehat  = matrix(ehat,ncol=sqrt(length(ehat)))
  Yhat  = muhat + ehat
  
  # Reflect monitored data
  Yhat2 <- Yhat
  Yhat2[seq(1,31,by=2),seq(1,31,by=2)]<-Ytrue[seq(1,31,by=2),seq(1,31,by=2)]
  
  # Boundary embedding
  Yhat3 <- Y
  Yhat3[3:nrow(Y),3:ncol(Y)] <- Yhat2
  
  out1 <- mean((Yhat3-Y)^2);if(out1 < 2) out[1] <- out1
  out[2] <- mean(evar)
  
  return(out)
}
# CalcMSPE.BSR(map,data0,NSBSR.PredMat); CalcMSPE.BSR(map,data0,PBSR.PredMat)

CalcMSPE.UK = function(map,geodata0,geodata00,alpha){
  
  out = c(NA,NA)
  
  Y = matrix(data0$Y,ncol=sqrt(map$nloc))
  Ytrue = Y[3:nrow(Y),3:ncol(Y)]
  
  fit0 = likfit(geodata00,ini.cov.pars=c(1,10),
                trend=trend.spatial(trend=~1+geodata00$covariate$X1+geodata00$covariate$X2,geodata=geodata00),
                cov.model="matern",nugget=0,fix.nugget=TRUE,messages=FALSE,kappa=alpha,fix.kappa=TRUE)
  krige = krige.control(type.krige = "ok",
                        trend.d = trend.spatial(trend=~1+geodata00$covariate$X1+geodata00$covariate$X2,geodata=geodata00),
                        trend.l = trend.spatial(trend=~1+geodata0$covariate$X1+geodata0$covariate$X2,geodata=geodata0),
                        cov.model = "matern",nugget=0,cov.pars=fit0$cov.pars,kappa=alpha)
  krige.UK = krige.conv(geodata=geodata00,locations=geodata0$coords,krige=krige)
  X = cbind(data0$X0,data0$X1,data0$X2)
  Ytilde = krige.UK$predict
  Ytilde = matrix(Ytilde,ncol=sqrt(nrow(data0)))

  out[1] <- mean((Y-Ytilde)^2)
  out[2] <- mean(krige.UK$krige.var)
  
  return(out)
}
# CalcMSPE.UK(map,geodata0,geodata00,0.1)

