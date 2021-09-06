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

# Extend.lambda = function(mat.lambda){
#   clambda = cbind(mat.lambda,mat.lambda[,1])
#   lambda0 = rbind(clambda,clambda[1,])
#   return(lambda0)
# }
# 
# Expand.iso = function(mat){
#   a = ncol(mat)-1
#   rmat = rbind(mat[1:a,],apply(mat,2,function(w) rev(w)))
#   cmat = cbind(rmat[,1:a],t(apply(t(rmat),2,function(w) rev(w))))
#   return(cmat)
# }

# Subfunctions for rho-grid search
Supp_grid_search_rho = function(simdata,rho.vec1,rho.vec2){
  n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
  freq.dist1 = d.o.o.w.x[1:mbar,1:mbar]
  freq.dist2 = d.o.o.w.y[1:mbar,1:mbar]
  inv.tau.dat =
    sapply(rho.vec1,function(rho1){
      sapply(rho.vec2,function(rho2){
        tau.mm2 = exp(-rho1*freq.dist1-rho2*freq.dist2)
        inv.tau.mm2 = solve(tau.mm2)
        return(inv.tau.mm2)
      })
    })
  out = matrix(inv.tau.dat,ncol=(length(rho.vec1)*length(rho.vec2)))
  return(out)
}

Effsolve = function(esti_cor_block,n){
  
  m = ncol(esti_cor_block[[1]])
  Zn = matrix(0,m,m)
  In = diag(rep(1,m))
  
  # P-matrix
  w = complex(real = cos(2*pi/n), imaginary = sin(2*pi/n))
  E = list()
  for(k in 1:n){
    wk = w^(k-1)
    E[[k]] <- In
    for(i in 2:n){
      E[[k]] <- rbind(E[[k]],wk^(i-1)*In)
    }
    if(k==1){P <- E[[k]]
    }else{P <- cbind(P,E[[k]])}
  }
  
  # S-matrix
  S = list()
  S.inv = list()
  for(q in 1:m){
    s = matrix(0,m,m)
    for(k in 1:n){
      s = s + w^((q-1)*(k-1))*(esti_cor_block[[k]]%>%as.matrix)
    }
    S[[q]] <- s
    S.inv[[q]] <- solve(s)
  }
  S.inv[[(n+1)]] <- Zn
  indmtx2 = matrix((n+1),n,n)-diag(rep((n+1),n))+diag(c(1:n))
  
  return( P%*%BlockToBCCB(S.inv,indmtx2)%*%Conj(P)/n )
  
}


# Effsolve = function(esti_cor_block,n){
#   
#   Zn = matrix(0,n,n)
#   In = diag(rep(1,n))
#   
#   # P-matrix
#   w = complex(real = cos(2*pi/n), imaginary = sin(2*pi/n))
#   E = list()
#   for(k in 1:n){
#     wk = w^(k-1)
#     E[[k]] <- In
#     for(i in 2:n){
#       E[[k]] <- rbind(E[[k]],wk^(i-1)*In)
#     }
#     if(k==1){P <- E[[k]]
#     }else{P <- cbind(P,E[[k]])}
#   }
#   
#   # S-matrix
#   S = list()
#   S.inv = list()
#   for(q in 1:n){
#     s = matrix(0,n,n)
#     for(k in 1:n){
#       s = s + w^((q-1)*(k-1))*(esti_cor_block[[k]]%>%as.matrix)
#     }
#     S[[q]] <- s
#     S.inv[[q]] <- solve(s)
#   }
#   S.inv[[(n+1)]] <- Zn
#   indmtx2 = matrix((n+1),n,n)-diag(rep((n+1),n))+diag(c(1:n))
#   
#   return( P%*%BlockToBCCB(S.inv,indmtx2)%*%Conj(P)/n )
#   
# }

iterate_spec <- function(y, observed, X = NULL, embed_fac = 1.2, burn_iters = 100, par_spec_fun = spec_AR1,
                         kern_parm = 2*pi/sqrt(sum(observed)), precond_method = "fft", m = 10,
                         silent = TRUE, max_iter = 200, tol = 1e-6,
                         converge_tol = 0.05,
                         ncondsim = 1){
  
  if(identical(par_spec_fun,FALSE)){
    do_parametric_filter <- FALSE
  } else {
    do_parametric_filter <- TRUE
  }
  
  # error for bad embedding factor value
  if( embed_fac < 1 ){
    stop("Embedding factor embed_fac cannot be less than 1")
  }
  
  # create embedded lattice
  nvec_obs <- dim(y)
  nvec <- round( nvec_obs*embed_fac )
  y_embed <- array(NA,nvec)
  y_embed[1:nvec_obs[1],1:nvec_obs[2]] <- y
  observed_embed <- array(FALSE,nvec)
  observed_embed[1:nvec_obs[1],1:nvec_obs[2]] <- observed
  
  # get setup for vecchia's approximation
  locsfull <- expand.grid( 1:nvec[1], 1:nvec[2] )
  locs <- locsfull[observed_embed,]
  NNarray <- findOrderedNN_kdtree(locs,m)
  NNarray[m+1,] <- (m+1):1 # need to make sure first part goes in the right order
  
  # get grid size and define kernel
  n <- prod(nvec)
  kern <- sqexp_kern(kern_parm, nvec)
  
  # do least squares to get initial estimate of mean
  if( !is.null(X) ){
    if( is.na(dim(X)[3]) ){
      n_covar <- 1
      dim(X) <- c(dim(X),1)
    } else {
      n_covar <- dim(X)[3]
    }
    X_embed <- array(NA, c(nvec,n_covar))
    for(j in 1:n_covar){
      X_embed[1:nvec_obs[1], 1:nvec_obs[2],j] <- X[ , , j]
    }
    X_embed_vecs <- X_embed
    dim(X_embed_vecs) <- c( prod(dim(X_embed)[1:2]), n_covar )
    infomat <- crossprod( X_embed_vecs[observed_embed,] )
    xcrossy <- crossprod( X_embed_vecs[observed_embed,], y_embed[observed_embed] )
    betahat <- solve( infomat, xcrossy )
    muvec <- X_embed_vecs %*% betahat
    mumat <- array( muvec, nvec )
    y0 <- y_embed - mumat
  } else {
    y0 <- y_embed - mean( y_embed, na.rm = TRUE )
    mumat <- array( mean( y_embed, na.rm = TRUE ), nvec )
    betahat <- NA
    infomat <- NA
  }
  
  # impute with mean
  y0[!observed_embed] <- mean(y0, na.rm = TRUE)
  
  # define likelihood function (written for AR1 filter)
  likfun <- function(x){
    expitx <- expit(x)/4
    return(-whittle_lik(pgram,expitx,par_spec_fun))
  }
  
  # set an initial value for optimization
  parm <- 1/8
  logitpar <- logit(parm*4)
  
  # loop over number of iterations
  
  for(k in 1:burn_iters){
    
    # periodogram
    pgram <- 1/n*abs( fft(y0) )^2
    # do the optimization
    if(do_parametric_filter){
      res <- optim(logitpar,likfun,method="Brent",lower=-6,upper=6)
      logitpar <- res$par         # update logit parameter
      parm <- expit(logitpar)/4   # update parameter
      # parametric approximation
      param_spec <- par_spec_fun(parm,nvec)
    } else {
      param_spec <- array(1,nvec)
    }
    # smooth ratio of periodogram to parametric approximation
    sm_pgram <- smooth_pgram(pgram/param_spec,kern, smoothlog = FALSE)
    # estimate is product of smoothed ratio and parametric spec
    spec <- sm_pgram*param_spec
    
    # every 10th iteration, re-estimate the mean
    if( k %% 10 == 0 ){
      if(!is.null(X)){
        Xsolve <- array(NA, c(nvec,n_covar))
        for(j in 1:n_covar){
          Xsolve[,,j] <- pcg_spec( X_embed[,,j], spec, observed_embed, precond_method, NNarray = NNarray, silent = silent )$x
        }
        Xsolve_vec <- Xsolve
        dim(Xsolve_vec) <- dim(X_embed_vecs)
        infomat <- crossprod( X_embed_vecs[observed_embed,], Xsolve_vec[observed_embed,] )
        xcrossy <- crossprod( Xsolve_vec[observed_embed,], y_embed[observed_embed] )
        betahat <- solve( infomat, xcrossy )
        muvec <- X_embed_vecs %*% betahat
        mumat <- array( muvec, nvec )
        y0 <- y_embed - mumat
      }
    }
    
    # do a conditional simulation with the new spectrum
    y0 <- condsim_spec(y = y0, spec = spec, obs = observed_embed, silent = silent, maxit = 500, precondmethod = precond_method, NNarray = NNarray, tol = tol)
  }
  
  spec_old <- spec
  
  for(k in 1:max_iter){
    pgram <- 1/n*abs( fft(y0) )^2
    # do the optimization
    if(do_parametric_filter){
      res <- optim(logitpar,likfun,method="Brent",lower=-6,upper=6)
      logitpar <- res$par         # update logit parameter
      parm <- expit(logitpar)/4   # update parameter
      # parametric approximation
      param_spec <- par_spec_fun(parm,nvec)
    } else {
      param_spec <- array(1,nvec)
    }
    # smooth ratio of periodogram to parametric approximation
    sm_pgram <- smooth_pgram(pgram/param_spec,kern, smoothlog = FALSE)
    # estimate is product of smoothed ratio and parametric spec
    spec <- sm_pgram*param_spec
    # update the spectrum estimate
    spec_new <- (k-1)/k*spec_old + 1/k*spec
    # compare difference to tolerance
    spec_sd <- sqrt( smooth_pgram( spec_old^2, kern^2 ) )
    
    criterion <- max( abs(spec_new - spec_old)/spec_sd )
    #if( k %% 10 == 0 ) cat(paste("Averaging Iteration",k,"Criterion =",round(criterion,4),"\n"))
    if( criterion < converge_tol ){
      break
    }
    
    # every 10th iteration, re-estimate the mean
    if( k %% 10 == 0 ){
      if(!is.null(X)){
        Xsolve <- array(NA, c(nvec,n_covar))
        for(j in 1:n_covar){
          Xsolve[,,j] <- pcg_spec( X_embed[,,j], spec_new, observed_embed, precond_method, NNarray = NNarray, silent = silent )$x
        }
        Xsolve_vec <- Xsolve
        dim(Xsolve_vec) <- dim(X_embed_vecs)
        infomat <- crossprod( X_embed_vecs[observed_embed,], Xsolve_vec[observed_embed,] )
        xcrossy <- crossprod( Xsolve_vec[observed_embed,], y_embed[observed_embed] )
        betahat <- solve( infomat, xcrossy )
        muvec <- X_embed_vecs %*% betahat
        mumat <- array( muvec, nvec )
        y0 <- y_embed - mumat
      }
    }
    
    spec_old <- spec_new
    y0 <- condsim_spec(y = y0, spec = spec_new, obs = observed_embed, silent = silent, maxit = 500, precondmethod = precond_method, NNarray = NNarray, tol = tol)
    
  }
  avg_spec <- spec_new
  
  # get estimate of the likelihood
  #cat("Computing estimate of likelihood \n")
  locsfull <- as.matrix( expand.grid( 1:nvec[1], 1:nvec[2] ) )
  locs <- locsfull[observed_embed,]
  NNarray <- findOrderedNN_kdtree(locs,40)
  covarray <- 1/prod(dim(avg_spec))*Re( fft( avg_spec, inverse = TRUE ) )
  yvec <- y_embed[observed_embed]
  loglik <- vecchiaLik(covarray,yvec,locs,NNarray)
  
  # get conditional expectation
  #cat("Computing Conditional Simulations \n")
  condexp <- condexp_spec(y0,avg_spec,observed_embed, silent=silent,
                          maxit=500,precondmethod = precond_method, NNarray = NNarray) +
    mumat
  
  # get conditional simulations
  condsim_array <- array( NA, c(nvec_obs,ncondsim) )
  #condsim_array <- array( NA, c(nvec,ncondsim) )
  for(j in 1:ncondsim){
    cursim <- condsim_spec(y0,avg_spec,observed_embed,silent=silent,
                           precondmethod=precond_method,NNarray=NNarray) +
      mumat
    condsim_array[,,j] <- cursim[1:nvec_obs[1],1:nvec_obs[2]]
    #condsim_array[,,j] <- cursim
  }
  
  return(list( spec = avg_spec, cov = covarray, loglik = loglik,
               condexp = condexp[1:nvec_obs[1],1:nvec_obs[2]], condsim = condsim_array,
               betahat = betahat, mumat = mumat[1:nvec_obs[1],1:nvec_obs[2]], betacov = solve(infomat)) )
}

spec_AR1 <- function(parm,nvec){
  
  n <- prod(nvec)
  filter0 <- array(0,nvec)
  filter0[1,1] = 1
  filter0[1,2] = -parm
  filter0[2,1] = -parm
  filter0[nvec[1],1] = -parm
  filter0[1,nvec[2]] = -parm
  
  fftfilter <- fft(filter0)
  spec <- 1/(abs(fftfilter)^2)
  spec <- n*spec/sum(spec)
  return(spec)
  
}

findOrderedNN_kdtree <- function(locs,m,mult=2,printsearch=FALSE){
  
  # number of locations
  n <- nrow(locs)
  
  # to store the nearest neighbor indices
  NNarray <- matrix(NA,n,m+1)
  # to the first mult*m+1 by brutce force
  NNarray[1:(mult*m+1),] <- findOrderedNN(locs[1:(mult*m+1),],m)
  
  query_inds <- (mult*m+2):n
  data_inds <- 1:n
  
  msearch <- m
  
  while( length(query_inds) > 0 ){
    
    msearch <- min( max(query_inds), 2*msearch )
    data_inds <- 1:max(query_inds)
    NN <- FNN::get.knnx( locs[data_inds,,drop=FALSE], locs[query_inds,,drop=FALSE], msearch )$nn.index
    less_than_k <- t(sapply( 1:nrow(NN), function(k) NN[k,] <= query_inds[k]  ))
    sum_less_than_k <- apply(less_than_k,1,sum)
    ind_less_than_k <- which(sum_less_than_k >= m+1)
    NN_less_than_k <- NN[ind_less_than_k,]
    
    NN_m <- t(sapply(ind_less_than_k,function(k) NN[k,][less_than_k[k,]][1:(m+1)] ))
    
    NNarray[ query_inds[ind_less_than_k], ] <- NN_m
    
    query_inds <- query_inds[-ind_less_than_k]
    if( printsearch ) print(length(query_inds))
    
  }
  
  return(NNarray)
}

findOrderedNN <- function( locs, m ){
  # find the m+1 nearest neighbors to locs[j,] in locs[1:j,]
  # by convention, this includes locs[j,], which is distance 0
  n <- dim(locs)[1]
  NNarray <- matrix(NA,n,m+1)
  for(j in 1:n ){
    distvec <- c(fields::rdist(locs[1:j,,drop=FALSE],locs[j,,drop=FALSE]) )
    NNarray[j,1:min(m+1,j)] <- order(distvec)[1:min(m+1,j)]
  }
  return(NNarray)
}

