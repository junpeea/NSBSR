
rm(list=ls())
setwd("~/Desktop/Code3(Acting)")
load("200601StudyDataset.RData")

X = as.data.frame(X.pred); colnames(X) <- c("1","R","S","NO2","VOC","UV","SO2");head(X)
Y = matrix(mydata$O3L3,ncol=1); colnames(Y) <- "O3L3"
simdata = list(D=grid.size, d=obs.delta, N=ND, nd=nd, nloc=nloc, nobs=nobs, map=map.grid,
               X=X, Y=Y)
mysimdata = data.frame(Y,X,simdata$map)

par(mfrow=c(1,2))
fields::image.plot(matrix(simdata$Y,(xrange+1),(yrange+1)))
# fields::image.plot(matrix(simdata$Y,(xrange+1),(yrange+1))[seq(1,xrange+1,by=2),seq(1,yrange+1,by=2)])
fields::image.plot(matrix(as.matrix(simdata$X) %*% fit1$coefficients,(xrange+1),(yrange+1)))

#######################################################################################################################
# Predictive Process (spLM)
#######################################################################################################################

library(spBayes)
library(dplyr)
library(mcmc)
library(coda)

time.pp <- system.time({
  
  
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  # fit0 = lm(O3L3~R+S+NO2+VOC+UV+SOl+SO2,data=train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  coords <- train.sim_data   %>% dplyr::select("x","y") %>% as.matrix
  # X <- train.sim_data %>% dplyr::select(R,S,NO2,VOC,UV,SOl,SO2) %>% cbind(1) %>% as.matrix
  X <- train.sim_data %>% dplyr::select(R,S,NO2,VOC,UV,SO2) %>% cbind(1) %>% as.matrix
  B <- as.matrix(fit0$coefficients)
  p <- length(B)
  
  sigma.sq <- 2
  sigma.sq <- 0.2
  sigma.sq <- 0.02
  tau.sq <- 0.1
  phi <- 3/0.5
  
  D <- as.matrix(dist(coords))
  R <- exp(-phi*D)
  w <- rmvn(1, rep(0,nrow(coords)), sigma.sq*R)
  y <- rnorm(nrow(coords), X%*%B + w, sqrt(tau.sq))
  
  n.samples <- 2500
  
  starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
  
  tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
  
  priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                   "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 0.1))
  
  priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
                   "sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))
  
  cov.model <- "exponential"
  
  n.report <- 100
  verbose <- TRUE
  
  m.1 <- spLM(y~X-1, coords=coords, starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=n.samples, verbose=verbose, n.report=n.report)
  
  # m.2 <- spLM(y~X-1, coords=coords, starting=starting,
  # tuning=tuning, priors=priors.2, cov.model=cov.model,
  # n.samples=n.samples, verbose=verbose, n.report=n.report)
  
  burn.in <- 0.5*n.samples
  
  X.ho = X.pred[obs.ind.mat2,]
  coords.pred = cbind(mydata$x,mydata$y)
  coords.ho = coords.pred[obs.ind.mat2,]
  
  m.1.pred <- spPredict(m.1, pred.covars=X.ho, pred.coords=coords.ho,start=0.5*n.samples)
  # m.1.pred <- spPredict(m.1, pred.covars=X.pred, pred.coords=coords.pred,start=0.5*n.samples)
  # m.2.pred <- spPredict(m.2, pred.covars=X.pred, pred.coords=coords.pred,start=0.5*n.samples)
  
  m.1.pred.summary <- summary(mcmc(t(m.1.pred[[1]])))$quantiles[,c(3,1,5)]
  # m.2.pred.summary <- summary(mcmc(t(m.2.pred[[1]])))$quantiles[,c(3,1,5)]
  
  PP.mat  = matrix(NA,(xrange+1),(yrange+1))
  PP.mat[obs.ind.mat2]  = m.1.pred.summary[,1]/sd(m.1.pred.summary[,1]) + X.ho %*% matrix(fit0$coefficients,ncol=1)
  int.mat = m.1.pred.summary[,3] - m.1.pred.summary[,2]
  
})

par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(PP.mat,zlim=c(270,350))

PP.mat[is.na(PP.mat)] = mean(PP.mat,na.rm=T)
Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(0,(xrange+1),(yrange+1))
sdhatmat[obs.ind.mat2] = (m.1.pred.summary[,3] - m.1.pred.summary[,2])/2
mean((Realdata >= PP.mat - 1.96 * sdhatmat) * (Realdata <= PP.mat + 1.96 * sdhatmat))


mae  = mean(abs((as.vector(O3L3.full) - as.vector(PP.mat))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(PP.mat) )^2,na.rm=T));round(rmse,3)
round(mean(as.vector(int.mat)/sd(m.1.pred.summary[,1])),3)
time.pp

setwd("~/Desktop/Code3(Acting)")
save.image(paste0(Sys.Date(),"_PP.Rdata"))

