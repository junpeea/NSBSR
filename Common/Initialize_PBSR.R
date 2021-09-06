########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : Initialize_PBSR.R

########################################################################################################
########################################################################################################


Initialize_PBSR = function(mcmcobj,simdata,kappa=0.5,nugget=0,range=10){
  
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  # data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  data00 = data.frame(simdata$map,simdata$X,simdata$Y) %>% filter(obs==4)
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  seednum <- mcmcobj$seed
  
  # initial beta
  XX = as.matrix(t(Xobs)%*%Xobs)
  betahat <- chol2inv(chol(XX)) %*% (t(Xobs) %*% Yobs)
  
  # perturbation
  seednum = seednum + 10; set.seed(seednum)
  betahat <- betahat + rnorm(length(betahat),0,5e-01)
  
  # initial residual
  resi <- Yobs - Xobs %*% betahat
  
  # initial sigma2
  sigma2 <- as.numeric(var(resi))
  
  # Parset : Parameter-set
  Parset = list()
  Parset$beta    <- betahat
  Parset$kappa   <- kappa
  Parset$sigma2  <- sigma2
  Parset$range   <- range
  Parset$nu2     <- nugget
  
  # MCMC storage
  n.beta = nrow(betahat);
  PBSR.EstiMat = as.mcmc(matrix(0,(mcmcobj$n.iter-mcmcobj$n.burn),(n.beta+1+1+1+1)))
  PBSR.PredMat = as.mcmc(matrix(0,(mcmcobj$n.iter-mcmcobj$n.burn),(simdata$nloc*3)))
  colnames(PBSR.EstiMat) <- c(paste('beta',c(0:(n.beta-1)),sep=''),'kappa','nugget','psill','range')
  colnames(PBSR.PredMat) <- c(paste('pred.mean',c(1:simdata$nloc),sep=''),paste('pred.error',c(1:simdata$nloc),sep=''),
                              paste('pred.vari',c(1:simdata$nloc),sep=''))
  
  out = list()
  out$n.chain <- mcmcobj$n.chain
  out$n.iter  <- mcmcobj$n.iter
  out$n.burn  <- mcmcobj$n.burn
  out$seed    <- seednum
  out$EstiMat <- PBSR.EstiMat
  out$PredMat <- PBSR.PredMat
  out$Parset  =  Parset
  
  return(out)
  
}

