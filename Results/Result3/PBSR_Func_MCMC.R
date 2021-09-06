########################################################################################################
########################################################################################################
# Thesis : Bayesian Spatial Process Regression with Nonparametric modeling of autocorrelation
# Project name : PBSR_func_MCMC

# Date : 201805xx
# Author : YB JUN / Chae Young Lim
########################################################################################################
########################################################################################################

########################################################################################################
### CONTENTS

### Initialize_PBSR(mcmcobj,simdata,kappa=0.5,nugget=0,range=10)
### Gibbs_PBSR(mcmcobj,simdata,beta.mu0=0,beta.sig20=100,a=1,b=10)
### Predict_PBSR = function(mcmcobj,simdata)
############################################################################################## CONTENTS
########################################################################################################


Initialize_PBSR = function(mcmcobj,simdata,kappa=0.5,nugget=0,range=10){
  
  data00  = simdata$DB[which(simdata$DB$obs==1),]
  Esti.grid = data00[,c("x","y","w.x","w.y")]
  X <- cbind(data00$X0,data00$X1,data00$X2); Y <- data00$Y
  seednum <- mcmcobj$seed
  
  # initial beta
  XX = as.matrix(t(X)%*%X)
  betahat <- chol2inv(chol(XX)) %*% (t(X) %*% Y)
  
  # perturbation
  seednum = seednum + 10; set.seed(seednum)
  betahat <- betahat + rnorm(length(betahat),0,5e-01)
  
  # initial residual
  resi <- Y - X %*% betahat
  
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
  PBSR.EstiMat = as.mcmc(matrix(0,(mcmcobj$n.iter-mcmcobj$n.burn),(simdata$n.beta+1+1+1+1)))
  PBSR.PredMat = as.mcmc(matrix(0,(mcmcobj$n.iter-mcmcobj$n.burn),(simdata$map$npred+simdata$map$npred+simdata$map$npred)))
  colnames(PBSR.EstiMat) <- c(paste('beta',c(0:(simdata$n.beta-1)),sep=''),'kappa','nugget','psill','range')
  colnames(PBSR.PredMat) <- c(paste('pred.mean',c(1:simdata$map$npred),sep=''),paste('pred.error',c(1:simdata$map$npred),sep=''),
                               paste('pred.vari',c(1:simdata$map$npred),sep=''))
  
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


Gibbs_PBSR = function(mcmcobj,simdata,beta.mu0=0,beta.sig20=100,a=100,b=100,fix.kappa=TRUE){
  
  new.mcmcobj <- mcmcobj
  
  data00  = simdata$DB[which(simdata$DB$obs==1),]
  geodata00 = as.geodata(data00,coords.col=1:2,data.col=13,covar.col=9:11)
  Esti.grid = data00[,c("x","y","w.x","w.y")]
  X <- cbind(data00$X0,data00$X1,data00$X2); Y <- data00$Y
  seednum <- mcmcobj$seed
  
  betahat <- mcmcobj$Parset$beta
  kappa   <- mcmcobj$Parset$kappa  # smoothness parameter
  sigma2  <- mcmcobj$Parset$sigma2 # partial sill
  range   <- mcmcobj$Parset$range  # range
  nu2     <- mcmcobj$Parset$nu2    # nugget
  
  resi    <- Y - X%*%betahat
  
  # Update beta
  temp = Esti.grid[,c("x","y")]
  Rstar = matrix(c(1,0,0,1),2,2)
  # Rstar = matrix(c(-1/sqrt(5),1/sqrt(5),0,1),2,2) %*% matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),2,2)
  # Rstar = matrix(c(-1/sqrt(101),1/sqrt(101),0,1),2,2) %*% matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),2,2)
  an.loc  = t(Rstar%*%t(temp))
  distvec = as.vector(as.matrix(dist(an.loc)))
  cor_n = matrix(Matern(distvec,smoothness=kappa,range=range),ncol=nrow(Esti.grid))
  cov_n = as.numeric(sigma2) * cor_n
  inv_cov_n = solve(cov_n)
  
  beta.sig2.star <- solve((t(X)%*%inv_cov_n%*%X + 1/beta.sig20))
  beta.mu.star   <- beta.sig2.star%*%(t(X)%*%inv_cov_n%*%Y + 1/beta.sig20*beta.mu0)
  seednum <- seednum + 10; set.seed(seednum)
  new.betahat <- as.matrix(mvrnorm(1,beta.mu.star,beta.sig2.star),length(beta.mu.star),1)
  
  # Update resi
  new.resi <- Y-X%*%new.betahat
  
  # Update sigma2
  tau  <- 1 / sigma2
  new.a = a + length(new.resi)/2
  new.b = b + 0.5*t(Re(new.resi))%*%inv_cov_n%*%Re(new.resi) 
  seednum <- seednum + 10; set.seed(seednum)
  new.tau = rgamma(1,shape=new.a,rate=new.b)
  new.sig2 = 1 / new.tau
  
  # Update range
  range.len = 101; range.len = 11
  range.vec = seq(1,simdata$map$gridsize,length.out=range.len)
  logprob0 = vector(len=(range.len-1),mode="numeric")
  for(k in 2:range.len){
    temp.cor_n = matrix(Matern(distvec,smoothness=kappa,range=range.vec[k]),ncol=nrow(Esti.grid))
    temp.cov_n = as.numeric(new.sig2) * temp.cor_n
    logprob0[k-1] = -log(range.vec[k]) + dmvnorm(t(new.resi),sigma=temp.cov_n,log=T)
    # cat(k-1, date(), "\n")
  }
  logprob = logprob0-max(logprob0)
  prob.v   = exp(logprob)/sum(exp(logprob))
  cumprob  = cumsum(prob.v[1:range.len])
  seednum <- seednum + 10; set.seed(seednum)
  uni <- runif(1)
  new.range  <- range.vec[which(cumprob>uni)][1]
  
  # Update nu2 (nugget assumed to be fixed zero)
  new.nu2 <- 0
  
  # Update kappa (fixed kappa version)
  new.kappa <- kappa
  
  # Update kappa
  if(fix.kappa==FALSE){
    kappa.vec = seq(0.1,2.0,by=0.1); kappa.len = length(kappa.vec)
    logprob0 = vector(len=kappa.len,mode="numeric")
    for(k in 1:kappa.len)
      logprob0[k] <- loglik.GRF(geodata = geodata00, obj.model = NULL, cov.model = "matern",
                               cov.pars=c(new.sig2,new.range), nugget = 0, kappa = kappa.vec[k],
                               trend=trend.spatial(trend=~1+geodata00$covariate$X1+geodata00$covariate$X2,geodata=geodata00))
    logprob = logprob0-max(logprob0)
    prob.v   = exp(logprob)/sum(exp(logprob))
    cumprob  = cumsum(prob.v[1:kappa.len])
    seednum <- seednum + 10; set.seed(seednum)
    uni <- runif(1)
    new.kappa  <- kappa.vec[which(cumprob>uni)][1]
  }
  
  new.mcmcobj$Parset$beta   <- new.betahat
  new.mcmcobj$Parset$kappa  <- new.kappa
  new.mcmcobj$Parset$sigma2 <- new.sig2
  new.mcmcobj$Parset$range  <- new.range
  new.mcmcobj$Parset$nu2    <- new.nu2
  new.mcmcobj$seed   <- seednum
  
  return(new.mcmcobj)
}

Predict_PBSR = function(mcmcobj,simdata){
  
  output = list()
  
  data0  = simdata$DB[which(simdata$DB$bdd2==0),]
  data00 = simdata$DB[which(simdata$DB$obs==1),]
  Esti.grid = data00[,c("x","y","w.x","w.y")]
  Pred.grid = data0[, c("x","y","w.x","w.y")]
  X <- cbind(data00$X0,data00$X1,data00$X2); Y <- data00$Y
  seednum <- mcmcobj$seed
  
  betahat = mcmcobj$Parset$beta
  kappa   = mcmcobj$Parset$kappa
  sig2    = mcmcobj$Parset$sigma2
  range   = mcmcobj$Parset$range
  nu2     = mcmcobj$Parset$nu2
  resi    = Y -X%*%betahat
  
  temp = Pred.grid[,c("x","y")]
  Rstar = matrix(c(1,0,0,1),2,2)
  # Rstar = matrix(c(-1/sqrt(5),1/sqrt(5),0,1),2,2) %*% matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),2,2)
  # Rstar = matrix(c(-1/sqrt(101),1/sqrt(101),0,1),2,2) %*% matrix(c(cos(pi/4),sin(pi/4),-sin(pi/4),cos(pi/4)),2,2)
  an.loc  = t(Rstar%*%t(temp))
  distvec = as.vector(as.matrix(dist(an.loc)))
  predcormtx = matrix(Matern(distvec,smoothness=kappa,range=range),ncol=nrow(Pred.grid))
  esticormtx = predcormtx[which(data0$obs==1),which(data0$obs==1)]
  inv_cov_n = solve(esticormtx)
  H.mtx = predcormtx[,which(data0$obs==1)]
  
  # Prediction results
  X = cbind(data0$X0,data0$X1,data0$X2)
  output$result1 = X%*%betahat
  output$result2 = (1/sqrt(sig2)) * H.mtx %*% inv_cov_n %*% resi
  # output$result3 = diag( pred_cor_mtx  - sqrt(tau.e)*H.mtx %*% cor_mtx_inv %*% t(H.mtx) )
  output$result3 = diag( predcormtx - H.mtx %*% inv_cov_n %*% t(H.mtx) )
  return(output)
  
}



