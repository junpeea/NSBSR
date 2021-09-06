########################################################################################################
########################################################################################################
# Thesis : Bayesian Spatial Process Regression with Nonparametric modeling of autocorrelation
# Project name : NSBSR_func_MCMC_ver4

# Date : 201811xx
# Author : YB JUN / Chae Young Lim
########################################################################################################
########################################################################################################

Update_beta_sig_phi_NSBSR = function(mcmcobj,simdata,beta.mu0=0,beta.sig20=100,a=100,b=100){
  
  new.mcmcobj <- mcmcobj
  
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  seednum <- mcmcobj$seed

  betahat  <- mcmcobj$Parset$beta
  tau.e    <- mcmcobj$Parset$tau.e
  lambda   <- mcmcobj$Parset$lambda
  mat.lambda = matrix(lambda,ncol=n)
  
  corr = fftw2d(mat.lambda,inverse=1)
  corr = matrix(Re(corr),ncol=n)
  
  # Estimataed correlation mtx
  esti_cor_block   = apply(corr,1,function(vec) circulant.spam(vec))
  esti_n = length(esti_cor_block)
  for(j in 1:esti_n){
    
    esti_cor_block[[j]] <- as.matrix(esti_cor_block[[j]])
    
  }
  cor_mtx_inv = Re(Effsolve(esti_cor_block,n))
  # indmtx = as.matrix(circulant.spam(seq(1,n,1)))
  # cor_mtx = BlockToBCCB(esti_cor_block,indmtx)
  # cor_mtx_inv = solve(cor_mtx)
  
  
  
  # Update beta
  beta.sig2.star <- solve(t(Xobs)%*%Re(cor_mtx_inv)%*%Xobs + 1/beta.sig20);
  temp = beta.sig2.star;  beta.sig2.star[lower.tri(temp)] <- beta.sig2.star[upper.tri(temp)]# numerical adjustment
  beta.mu.star   <- beta.sig2.star%*%(t(Xobs)%*%Re(cor_mtx_inv)%*%Yobs + 1/beta.sig20*beta.mu0*rep(1,length(betahat)))
  seednum <- seednum + 10; set.seed(seednum)
  if(is.positive.definite(beta.sig2.star)){
    new.betahat <- matrix(mvrnorm(1,beta.mu.star,beta.sig2.star),length(beta.mu.star),ncol=1)
  }else{
    cat("positive definite failed","\n")
    new.betahat <- matrix(mvrnorm(1,beta.mu.star,solve(t(X)%*%X)),length(beta.mu.star),ncol=1)
  }
  
  # Update sig.e
  new.resi <- Yobs - Xobs%*%new.betahat
  sc.resi <- new.resi / sd(new.resi)
  new.a <- a + n2/2
  new.b <- b + 0.5*t(sc.resi)%*%cor_mtx_inv%*%sc.resi
  seednum <- seednum + 10; set.seed(seednum)
  new.tau.e <- rgamma(1,shape=new.a,rate=new.b)
  if(is.na(new.tau.e)) new.tau.e <- as.numeric(1/var(new.resi))
  
  # Update phi
  sc.resi    <- new.resi * sqrt(new.tau.e)
  new.phi    <- CompLogPerio_fft(sc.resi,mode=2)
  
  new.mcmcobj$Parset$beta   <- new.betahat
  new.mcmcobj$Parset$tau.e  <- new.tau.e
  new.mcmcobj$Parset$phi    <- new.phi
  new.mcmcobj$seed   <- seednum
  
  return(new.mcmcobj)
}


Update_Theta_NSBSR = function(mcmcobj,simdata){
  
  new.mcmcobj <- mcmcobj
  
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  seednum <- mcmcobj$seed
  
  phi     <- mcmcobj$Parset$phi
  nu.th   <- mcmcobj$Parset$nu.th
  tau.th  <- mcmcobj$Parset$tau.th
  rho.th1 <- mcmcobj$Parset$rho.th1
  rho.th2 <- mcmcobj$Parset$rho.th2
  label   <- mcmcobj$Parset$label
  klabel  <- mcmcobj$Parset$klabel
  vlabel  <- mcmcobj$Parset$vlabel
  
  # Direct evaluation -> storage limitation
  freq.dist1 = d.o.o.w.x[1:mbar,1:mbar]
  freq.dist2 = d.o.o.w.y[1:mbar,1:mbar]
  tau.star  = solve(tau.th*solve(exp(-rho.th1*freq.dist1-rho.th2*freq.dist2))+diag(1/vlabel))
  
  nu.star  = tau.star%*%diag(1/vlabel)%*%(phi-klabel-nu.th)+nu.th
  seednum <- seednum + 10; set.seed(seednum)
  theta = matrix(rmvnorm(1,mean=nu.star,sigma=tau.star),ncol=1)
  
  # lambda: rescaled to 1
  if(n%%2==0){lambda0 = exp(c(theta,rev(theta)))
  }else{lambda0 = exp(c(theta,rev(theta[1:(length(theta)-1)])))} 
  newlambda0 = lambda0 / sum(lambda0)
  newtheta = log(newlambda0)[1:mbar]
  
  # output
  new.mcmcobj$Parset$lambda <- newlambda0
  new.mcmcobj$Parset$theta  <- newtheta
  new.mcmcobj$seed   <- seednum
  
  return(new.mcmcobj)
}



Update_Theta_prior_NSBSR = function(mcmcobj,simdata,rho.list,c=1,d=1){
  
  new.mcmcobj <- mcmcobj
  
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  seednum <- mcmcobj$seed
  
  theta   <- mcmcobj$Parset$theta
  nu.th   <- mcmcobj$Parset$nu.th
  tau.th  <- mcmcobj$Parset$tau.th
  rho.th1  <- mcmcobj$Parset$rho.th1
  rho.th2  <- mcmcobj$Parset$rho.th2
  
  # Update tau.th
  freq.dist1 = d.o.o.w.x[1:mbar,1:mbar]
  freq.dist2 = d.o.o.w.y[1:mbar,1:mbar]
  tau.mm=(1/tau.th)*exp(-rho.th1*freq.dist1-rho.th1*freq.dist2)
  
  new.c <- c + m/2
  new.d <- d + 0.5*t(theta-nu.th)%*%chol2inv(chol(tau.mm))%*%(theta-nu.th)
  seednum <- seednum + 10; set.seed(seednum)
  new.tau.th = rgamma(1,shape=new.c,rate=new.d)
  if(is.na(new.tau.th)) new.tau.th <- 1
  
  # Update rho.th1, rho.th2
  rho.vec = rho.list[[1]]; rho.grid = rho.list[[2]]
  logprob0 = apply(rho.grid,2,function(vec){
    inv.tau.mm2 = new.tau.th*matrix(vec,nrow=mbar)
    logprob0 = -0.5*t(theta)%*%inv.tau.mm2%*%theta
    logprob0 = logprob0 - 0.5*determinant(inv.tau.mm2)$mod
    return(as.numeric(logprob0))
  })
  
  logprob  = logprob0-max(logprob0)
  prob.v   = exp(logprob)/sum(exp(logprob))
  cumprob  = cumsum(prob.v)
  seednum <- seednum + 10; set.seed(seednum)
  uni <- runif(1); r = min(which(cumprob>uni)) 
  # new.rho.th1 <- rho.vec[(r%%100)+1]; new.rho.th2 <- rho.vec[ceiling(r/100)]
  # new.rho.th1 <- rho.vec[1]; new.rho.th2 <- rho.vec[1]
  temp = expand.grid(c(1:length(rho.vec)),c(1:length(rho.vec)))
  new.rho.th1 <- rho.vec[temp[r,1]]; new.rho.th2 <- rho.vec[temp[r,2]]
  if(is.na(new.rho.th1) | is.na(new.rho.th2)){new.rho.th1 <- rho.vec[1];new.rho.th2 <- rho.vec[1];} 
  
  # output
  new.mcmcobj$Parset$tau.th  <- new.tau.th
  new.mcmcobj$Parset$rho.th1 <- new.rho.th1
  new.mcmcobj$Parset$rho.th2 <- new.rho.th2
  new.mcmcobj$seed   <- seednum
  
  return(new.mcmcobj)  
}


Update_Mixture_NSBSR = function(mcmcobj,simdata){
  
  kl.w0 = c(-4.63,-2.87,-1.44,-0.33,0.76)
  kl.wj = c(-2.20,-0.80,-0.55,-0.035,0.48)
  vl.w0 = c(8.75,1.95,0.88,0.45,0.41)
  vl.wj = c(1.93,1.01,0.69,0.60,0.29)
  pl.w0 = c(0.13,0.16,0.23,0.22,0.25)
  pl.wj = c(0.19,0.11,0.27,0.25,0.18)
  
  new.mcmcobj <- mcmcobj
  
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  seednum <- mcmcobj$seed
  
  phi    <- mcmcobj$Parset$phi
  theta  <- mcmcobj$Parset$theta
  xi = phi - theta
  label  = rep(0,m^2);klabel = rep(0,m^2);vlabel  = rep(0,m^2)
  
  bddind = simdata$map %>% filter(obs==4) %>% dplyr::select(bdd)
  kl.mtx = matrix(rep(kl.wj,n2),ncol=n2); kl.mtx[,which(bddind==1)] <- kl.w0
  vl.mtx = matrix(rep(vl.wj,n2),ncol=n2); vl.mtx[,which(bddind==1)] <- vl.w0
  pl.mtx = matrix(rep(pl.wj,n2),ncol=n2); pl.mtx[,which(bddind==1)] <- pl.w0
  
  for(i in 1:mbar){
    log.alpha0 = array(0,5)
    for(j in 1:5){
      seednum = seednum + 10; set.seed(seednum)
      log.alpha0[j] <- log(pl.mtx[j,i]) + dnorm(xi[i],mean=kl.mtx[j,i],sd=sqrt(vl.mtx[j,i]),log=TRUE)
    }
    log.alpha = log.alpha0 - max(log.alpha0)
    pi.hat = exp(log.alpha)/sum(exp(log.alpha))
    label[i]  <- sample(c(1:5),size=1,prob=pi.hat)
    klabel[i] <- kl.mtx[label[i],i]
    vlabel[i] <- vl.mtx[label[i],i]
  }
  
  # output
  new.mcmcobj$Parset$label  <- label
  new.mcmcobj$Parset$klabel <- klabel
  new.mcmcobj$Parset$vlabel <- vlabel
  new.mcmcobj$seed   <- seednum
  
  return(new.mcmcobj)
  
}


Predict_NSBSR = function(mcmcobj,Esti.summary,simdata){
  
  output = list()
  
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  data0  = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y)
  data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
  Esti.grid = data00 %>% dplyr::select("x","y","w.x","w.y")
  Pred.grid = data0  %>% dplyr::select("x","y","w.x","w.y")
  Xobs <- data00 %>% dplyr::select(X0,X1,X2) %>% as.matrix; Yobs <- data00$simdata.Y
  seednum <- mcmcobj$seed
  
  betahat  <- mcmcobj$Parset$beta
  tau.e    <- mcmcobj$Parset$tau.e
  lambda   <- mcmcobj$Parset$lambda
  resi    = Yobs -Xobs%*%betahat
  
  meantheta = Esti.summary$statistics[(mcmcobj$n.beta+1+1):(mcmcobj$n.beta+1+mcmcobj$n.theta),1]
  if(n%%2==0){
    mat.theta  = matrix(c(meantheta,rev(meantheta)),ncol=n)
  }else{
    mat.theta  = matrix(c(meantheta,rev(meantheta[1:(length(meantheta)-1)])),ncol=n)
  } 
  mat.lambda = exp(mat.theta)
  
  # estimated correlation
  corr = Re(fftw2d(mat.lambda,inverse=1))

  # corr = as.vector(corr)
  # cor_mtx = as.matrix(toeplitz.spam(corr))
  
  #############################################################################################################
  
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
  cor_mtx_inv = solve(cor_mtx)
  
  # temp = toeplitz.spam(c(1:17))%>%as.matrix
  # temp2 = lapply(c(1:17),function(w) toeplitz.spam(corr[,w]))
  # temp5 = kronecker(temp,matrix(1,17,17))
  # cor_mtx = matrix(0,289,289)
  # for(i in 1:17){
  #   
  #   temp6 = matrix(as.numeric(temp5 %in% i),nrow=289)
  #   temp7 = kronecker(matrix(1,17,17),temp2[[i]] %>% as.matrix)
  #   cor_mtx = cor_mtx + temp6*temp7
  #   
  # }
  # cor_mtx_inv = solve(cor_mtx)
  
  #############################################################################################################
  
  # correlation interpolation
  curr.obj <- list( x=seq(0,simdata$D[1],by=simdata$d[1]), y=seq(0,simdata$D[2],by=simdata$d[2]), z= corr )
  make.surface.grid( list(seq(0,simdata$D[1]),seq(0,simdata$D[2])) )-> loc
  interp.surface(curr.obj, loc)-> look
  pred_corr = matrix(look,ncol=sqrt(length(look)),byrow=T)
  
  
  temp = toeplitz.spam(c(1:33))%>%as.matrix
  temp2 = lapply(c(1:33),function(w) toeplitz.spam(pred_corr[,w]))
  pred_cor_mtx = matrix(NA,33^2,33^2)
  for(i in 1:(33^2)){
    for(j in 1:(33^2)){
      a = i%%33; b = ceiling(i/33)
      c = j%%33; d = ceiling(j/33)
      if(a==0) a=33; if(c==0) c=33;
      pred_cor_mtx[i,j] = temp2[[temp[a,c]]][b,d]
    }
  }
  
  # temp = toeplitz.spam(c(1:33))%>%as.matrix
  # temp2 = lapply(c(1:33),function(w) toeplitz.spam(pred_corr[,w]))
  # temp5 = kronecker(temp,matrix(1,33,33))
  # pred_cor_mtx = matrix(0,33^2,33^2)
  # for(i in 1:33){
  #   
  #   temp6 = matrix(as.numeric(temp5 %in% i),nrow=,33^2)
  #   temp7 = kronecker(matrix(1,33,33),temp2[[i]] %>% as.matrix)
  #   pred_cor_mtx = pred_cor_mtx + temp6*temp7
  #   
  # }
  
  # pred_corr = as.vector(pred_corr)
  # pred_cor_mtx = as.matrix(toeplitz.spam(pred_corr))
  
  H.mtx = pred_cor_mtx[,which(data0$obs==4)]
  
  # Prediction results
  X = data0 %>% dplyr::select(X0,X1,X2) %>% as.matrix
  output$result1 = X%*%betahat
  output$result2 = sqrt(tau.e)*H.mtx %*% cor_mtx_inv %*% resi
  output$result3 = diag( pred_cor_mtx  - H.mtx %*% cor_mtx_inv %*% t(H.mtx) )
  output$result4 = cor_mtx[1:17,1:17]
  return(output)
  
}


CompLogPerio_fft = function(sc.resi,mode=1){
  n2 = length(sc.resi); mbar=ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  e = matrix(sc.resi,ncol=n)
  D = fftw2d(e); D2 = Re(D*Conj(D))
  lambda <- as.vector(D2/sum(D2))
  temp <- matrix(log(lambda),ncol=n)
  if(mode==1) new.phi <- as.vector(temp[1:m,1:m])
  if(mode==2) new.phi <- as.vector(temp)[1:mbar]
  return(new.phi)
}


LambdaToGammanij = function(u1,u2,lambda){
  n2 = length(lambda)
  n = sqrt(n2) ; m = floor(n/2)
  w = complex(real = cos(2*pi/n),imaginary = sin(2*pi/n))
  s = 0 ; iter = 1
  K = as.vector(matrix(rep(seq(-m,by=1,length.out=n),n),n))
  J = as.vector(t(matrix(rep(seq(-m,by=1,length.out=n),n),n)))
  return(sum(lambda*w^(u1*K+u2*J))/n^2)
  # for(j in seq(-m,by=1,length.out=n)){
  #   for(k in seq(-m,by=1,length.out=n)){
  #     s = s + lambda[iter]*w^(u1*k+u2*j)
  #     iter = iter + 1
  #   }
  # }
  # return(s/n^2)
}

LambdaToGammaBlock = function(lambda){
  
  n2 = length(lambda); n = sqrt(n2)
  u1vec = seq(0,(n-1),1)
  
  blockmtx = list()
  for(u2 in 0:(n-1)){
    f1 = function(u1){
      return(LambdaToGammanij(u1,u2,lambda))
    }
    realpart <- as.matrix(circulant.spam(Re(sapply(u1vec,FUN = f1))))
    imagpart <- as.matrix(circulant.spam(Im(sapply(u1vec,FUN = f1))))
    blockmtx[[(u2+1)]] = matrix(complex(real = realpart,imaginary = imagpart),ncol=n)
  }
  
  return(blockmtx)
  
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

LambdaToGamman = function(lambda){
  
  n2 = length(lambda); n = sqrt(n2)
  u1vec = u1vec = seq(0,(n-1),1)
  indmtx = as.matrix(circulant.spam((u1vec+1)))
  blocklist = LambdaToGammaBlock(lambda)
  finalmtx = BlockToBCCB(blocklist,indmtx)
  return(finalmtx)
  
}


