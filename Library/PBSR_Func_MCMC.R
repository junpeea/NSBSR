########################################################################################################
########################################################################################################
# Thesis : Bayesian Spatial Process Regression with Nonparametric modeling of autocorrelation
# Project name : PBSR_func_MCMC

# Date : 201811xx
# Author : YB JUN / Chae Young Lim
########################################################################################################
########################################################################################################


Gibbs_PBSR = function(mcmcobj,simdata,beta.mu0=0,beta.sig20=100,a=100,b=100,fix.kappa=TRUE){
  
  new.mcmcobj <- mcmcobj
  
  # data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  data00 = data.frame(simdata$map,simdata$X,simdata$Y) %>% filter(obs==4)
  summary(simdata$Y)
  table(simdata$map$obs)
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  seednum <- mcmcobj$seed
  # geodata00 = as.geodata(data00,coords.col=1:2,data.col=11,covar.col=7:9)
  geodata00 = as.geodata(data00,coords.col=1:2,data.col=10,covar.col=7:9)
  
  betahat <- mcmcobj$Parset$beta
  kappa   <- mcmcobj$Parset$kappa  # smoothness parameter
  sigma2  <- mcmcobj$Parset$sigma2 # partial sill
  range   <- mcmcobj$Parset$range  # range
  nu2     <- mcmcobj$Parset$nu2    # nugget
  
  resi    <- Yobs - Xobs%*%betahat
  
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
  
  beta.sig2.star <- solve((t(Xobs)%*%inv_cov_n%*%Xobs + 1/beta.sig20))
  beta.mu.star   <- beta.sig2.star%*%(t(Xobs)%*%inv_cov_n%*%Yobs + 1/beta.sig20*beta.mu0)
  seednum <- seednum + 10; set.seed(seednum)
  new.betahat <- as.matrix(mvrnorm(1,beta.mu.star,beta.sig2.star),length(beta.mu.star),1)
  
  # Update resi
  new.resi <- Yobs-Xobs%*%new.betahat
  
  # Update sigma2
  tau  <- 1 / sigma2
  new.a = a + length(new.resi)/2
  new.b = b + 0.5*t(Re(new.resi))%*%inv_cov_n%*%Re(new.resi) 
  seednum <- seednum + 10; set.seed(seednum)
  new.tau = rgamma(1,shape=new.a,rate=new.b)
  new.sig2 = 1 / new.tau
  
  # Update range
  range.len = 101; range.len = 11
  range.vec = seq(1,simdata$D[1],length.out=range.len)
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
  
  # data0  = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y)
  data0  = data.frame(simdata$map,simdata$X,simdata$Y)
  # data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  data00 = data.frame(simdata$map,simdata$X,simdata$Y) %>% filter(obs==4)
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  Pred.grid = data0  %>% dplyr::select("x","y","w.x","w.y")
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  seednum <- mcmcobj$seed
 
  
  betahat = mcmcobj$Parset$beta
  kappa   = mcmcobj$Parset$kappa
  sig2    = mcmcobj$Parset$sigma2
  range   = mcmcobj$Parset$range
  nu2     = mcmcobj$Parset$nu2
  resi    = Yobs -Xobs%*%betahat
  
  temp = Pred.grid[,c("x","y")]
  Rstar = matrix(c(1,0,0,1),2,2)
  an.loc  = t(Rstar%*%t(temp))
  distvec = as.vector(as.matrix(dist(an.loc)))
  predcormtx = matrix(Matern(distvec,smoothness=kappa,range=range),ncol=nrow(Pred.grid))
  esticormtx = predcormtx[which(data0$obs==4),which(data0$obs==4)]
  inv_cov_n = solve(esticormtx)
  H.mtx = predcormtx[,which(data0$obs==4)]
  
  # Prediction results
  TX = cbind(data0$X0,data0$X1,data0$X2)
  output$result1 = TX%*%betahat
  output$result2 = (1/sqrt(sig2)) * H.mtx %*% inv_cov_n %*% resi
  output$result3 = diag( predcormtx - H.mtx %*% inv_cov_n %*% t(H.mtx) )
  output$result4 = esticormtx[1:64,1:64]
  return(output)
  
}
