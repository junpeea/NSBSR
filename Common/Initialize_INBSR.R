########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : Initialize_NSBSR.R

########################################################################################################
########################################################################################################

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


Initialize_NSBSR = function(mcmcobj,simdata){
  
  kl.w0 = c(-4.63,-2.87,-1.44,-0.33,0.76)
  kl.wj = c(-2.20,-0.80,-0.55,-0.035,0.48)
  vl.w0 = c(8.75,1.95,0.88,0.45,0.41)
  vl.wj = c(1.93,1.01,0.69,0.60,0.29)
  pl.w0 = c(0.13,0.16,0.23,0.22,0.25)
  pl.wj = c(0.19,0.11,0.27,0.25,0.18)
  
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  data00 = data.frame(simdata$map,simdata$X,simdata$Y) %>% filter(obs > 0)
  # Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  Xobs <- data00 %>% dplyr::select("X1","R","S","NO2","VOC","UV","SO2") %>% as.matrix; Yobs <- data00$simdata.Y
  Yobs[which(data00$obs==2)] <- NA # unobserved monitors
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  seednum <- mcmcobj$seed
  
  # initial Y
  Yobs[which(is.na(Yobs)==1)] <- mean(Yobs,na.rm=T) # simple imputation
  
  # initial beta
  XX = as.matrix(t(Xobs)%*%Xobs)
  betahat <- chol2inv(chol(XX)) %*% (t(Xobs) %*% Yobs)
  
  # perturbation
  seednum = seednum + 10; set.seed(seednum)
  betahat <- betahat + rnorm(length(betahat),0,5e-01)
  
  # initial residual (resi, sc.resi)
  resi    <- Yobs - Xobs %*% betahat
  sig2.e  <- as.numeric(var(resi))
  sc.resi <- resi / sqrt(sig2.e)
  
  # initial lambda
  e = matrix(sc.resi,ncol=n)
  # D = fftw2d(e); D2 = Re(D*Conj(D))
  D = t(mvfft(t(mvfft(e)))); D2 = Re(D*Conj(D))
  new.lambda <- as.vector(D2/sum(D2))
  
  # initial phi
  new.phi <- log(new.lambda)[1:mbar]
  
  # initial theta  ( kernel smoothing lambda)
  smooth.lambda <- ksmooth(x=c(1:n^2),y=new.lambda,n.points=n^2,bandwidth=3)$y
  new.smooth.lambda <- smooth.lambda / sum(smooth.lambda)
  new.theta = log(new.smooth.lambda)[1:mbar]
  
  
  # initial xi (label-par.)
  xi = new.phi - new.theta
  label  = rep(0,mbar);klabel = rep(0,mbar);vlabel = rep(0,mbar)
  kl.mtx = matrix(rep(kl.wj,n2),ncol=n2); kl.mtx[,which(data00$bdd==1)] <- kl.w0
  vl.mtx = matrix(rep(vl.wj,n2),ncol=n2); vl.mtx[,which(data00$bdd==1)] <- vl.w0
  pl.mtx = matrix(rep(pl.wj,n2),ncol=n2); pl.mtx[,which(data00$bdd==1)] <- pl.w0
  for(i in 1:mbar){
    log.alpha0 = array(0,5)
    for(j in 1:5){
      seednum = seednum + 1; set.seed(seednum)
      log.alpha0[j] <- log(pl.mtx[j,i]) + dnorm(xi[i],mean=kl.mtx[j,i],sd=sqrt(vl.mtx[j,i]),log=TRUE)
    }
    log.alpha = log.alpha0 - max(log.alpha0)
    pi.hat = exp(log.alpha)/sum(exp(log.alpha))
    label[i]  <- sample(c(1:5),size=1,prob=pi.hat)
    klabel[i] <- kl.mtx[label[i],i]
    vlabel[i] <- vl.mtx[label[i],i]
  }
  
  # initial nu
  nu.th = rep(0,mbar)
  # library(spectralGP)
  # n1=simdata$D[1]; n2=simdata$D[2]; alpha=1/simdata$d[1]
  # gp2=gp(c(n1,n2),matern.specdens,c(0.4,alpha))
  # dens2=matern.specdens(gp2$omega,c(0.4,alpha),d=2)
  # dens2mat = matrix(dens2,n1,n2)
  # logdens2mat = log(dens2mat / sum(dens2mat))
  # nu.th = as.vector(logdens2mat[seq(2,32,by=2),seq(2,32,by=2)])[1:mbar]
  # nu.th = as.vector(logdens2mat[seq(2,32,by=2),seq(2,32,by=2)])[1:mbar]
  # temp = matrix(c(nu.th,rev(nu.th)),16,16)
  # temp = matrix(0,16,16)
  # temp[5:12,] <- -10; temp[,5:12] <- -10
  # nu.th = as.vector(temp)[1:mbar]
  
  # temp = matrix(c(new.theta[1:(length(new.theta)-1)],rev(new.theta)),n,n)
  # y = as.vector(temp)
  # d = expand.grid(c(1:n),c(1:n))
  # temp2 = smooth.2d(Y=y,x=d,theta=2.5,nrow=(2*n),ncol=(2*n))
  # temp3 = temp2$z[seq(1,(2*n),by=2),seq(1,(2*n),by=2)]
  # temp4 = as.vector(temp3)
  # nu.th = temp4[1:mbar]*0.05
  
  
  # Parset : Parameter-set
  Parset = list()
  Parset$beta    <- betahat
  Parset$tau.e  <- (1/sig2.e)
  Parset$phi     <- new.phi
  Parset$lambda  <- as.vector(new.lambda)
  Parset$theta   <- new.theta
  Parset$nu.th   <- nu.th
  Parset$tau.th  <- 1
  Parset$rho.th1 <- 1
  Parset$rho.th2 <- 1
  Parset$label   <- label
  Parset$klabel  <- klabel
  Parset$vlabel  <- vlabel
  
  # MCMC storage
  n.beta = nrow(betahat); n.theta = length(new.theta)
  NSBSR.EstiMat = as.mcmc(matrix(NA,(mcmcobj$n.iter-mcmcobj$n.burn),(n.beta+1+n.theta+1+1)))
  NSBSR.PredMat = as.mcmc(matrix(0,(mcmcobj$n.iter-mcmcobj$n.burn),(simdata$nloc*3)))
  colnames(NSBSR.EstiMat) <- c(paste('beta',c(0:(n.beta-1)),sep=''),'tau.e',paste('theta',c(1:n.theta),sep=''),'rho1','rho2')
  colnames(NSBSR.PredMat) <- c(paste('pred.mean',c(1:simdata$nloc),sep=''),paste('pred.error',c(1:simdata$nloc),sep=''),
                               paste('pred.vari',c(1:simdata$nloc),sep=''))
  
  out = list()
  out$n.chain <- mcmcobj$n.chain
  out$n.iter  <- mcmcobj$n.iter
  out$n.burn  <- mcmcobj$n.burn
  out$n.beta  <- n.beta
  out$n.theta <- n.theta
  out$seed    <- seednum
  out$EstiMat <- NSBSR.EstiMat
  out$PredMat <- NSBSR.PredMat
  out$Parset  =  Parset
  out$sampleY <- Yobs
  
  return(out)
  
}


