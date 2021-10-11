
#######################################################################################################################
# Predictive Process (spLM)
#######################################################################################################################

library(spBayes)

time.pp <- system.time({
  
  
  rmvn <- function(n, mu=0, V = matrix(1)){
    p <- length(mu)
    if(any(is.na(match(dim(V),p))))
      stop("Dimension problem!")
    D <- chol(V)
    t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
  }
  
  set.seed(1)
  
  ### Generate process and data
  n <- 30
  temp = expand.grid(c(0:(n-1)),c(0:(n-1)))
  sim_process = data.frame(x=temp$Var1,y=temp$Var2)
  sim_process$Y = as.vector(Y0)
  
  scale = function(x) {(x-mean(x))/sd(x)}
  
  par(mfrow=c(1,3))
  fields::image.plot(Realdata)
  fields::image.plot(matrix(sim_process$Y,30,30))
  fields::image.plot(matrix(Y00,30,30))
  
  ind = is.na(as.vector(Y00));summary(ind)
  sim_data <- sim_process[which(ind==FALSE),]
  sim_data$x1 = rep(1,nrow(sim_data))
  sim_data$x2 = scale(sim_process$x[which(ind==FALSE)])
  sim_data$x3 = scale(sim_process$y[which(ind==FALSE)])
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  sim_data$z = fit0$residuals
  
  head(sim_data)
  coordinates(sim_data) = ~x + y # change into an sp object
  
  coords <- cbind(sim_data$x,sim_data$y)
  X <- cbind(sim_data$x1,sim_data$x2,sim_data$x3)
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
  # y <- rnorm(nrow(coords), X%*%B + w, sqrt(tau.sq))
  y <- rnorm(nrow(coords), sim_data$Y, sqrt(tau.sq))
  
  n.samples <- 10000
  
  starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
  
  tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
  
  priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
                   "phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
                   "tau.sq.IG"=c(2, 0.1))
  
  priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
                   "sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))
  
  cov.model <- "exponential"
  
  n.report <- 1000
  verbose <- TRUE
  
  m.1 <- spLM(y~X-1, coords=coords, starting=starting,
              tuning=tuning, priors=priors.1, cov.model=cov.model,
              n.samples=n.samples, verbose=verbose, n.report=n.report)
  
  m.2 <- spLM(y~X-1, coords=coords, starting=starting,
  tuning=tuning, priors=priors.2, cov.model=cov.model,
  n.samples=n.samples, verbose=verbose, n.report=n.report)
  
  burn.in <- 0.5*n.samples
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X.pred = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  coords.pred = cbind(sim_process$x,sim_process$y)
  m.1.pred <- spPredict(m.1, pred.covars=X.pred, pred.coords=coords.pred,start=0.5*n.samples)
  # m.2.pred <- spPredict(m.2, pred.covars=X.pred, pred.coords=coords.pred,start=0.5*n.samples)
  
  m.1.pred.summary <- summary(mcmc(t(m.1.pred[[1]])))$quantiles[,c(3,1,5)]
  # m.2.pred.summary <- summary(mcmc(t(m.2.pred[[1]])))$quantiles[,c(3,1,5)]
  
  PP.mat  = matrix(m.1.pred.summary[,1],30,30)
  fields::image.plot(PP.mat)
  int.mat = m.1.pred.summary[,3] - m.1.pred.summary[,2]
  mean((Realdata > matrix(m.1.pred.summary[,2],30,30)) * (Realdata < matrix(m.1.pred.summary[,3],30,30)))
  
  par(mfrow=c(1,2))
  fields::image.plot(Realdata)
  fields::image.plot(matrix(X.pred%*%fit0$coefficient,30,30))
  
  PP.mat  = matrix(m.2.pred.summary[,1],30,30)
  fields::image.plot(PP.mat)
  int.mat = m.2.pred.summary[,3] - m.2.pred.summary[,2]
  mean((Realdata > matrix(m.2.pred.summary[,2],30,30)) * (Realdata < matrix(m.2.pred.summary[,3],30,30)))
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  pred.NSBSR       = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  PP.mat  = matrix(m.1.pred.summary[,1],30,30)
  int.mat = m.1.pred.summary[,3] - m.1.pred.summary[,2]
  
  summary((t(mylist[[k]])[1:30,1:30] - t(PP.mat[1:30,1:30]) )^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(PP.mat[1:30,1:30])%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(PP.mat,zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(PP.mat)
  
  summary((t(mylist[[k]])[1:30,1:30] - t(pred.NSBSR[1:30,1:30]))^2 %>% as.vector)
  cor(x=t(mylist[[k]])[1:30,1:30]%>%as.vector,y=t(pred.NSBSR[1:30,1:30])%>%as.vector,use="complete.obs")
  # summary((t(mylist[[1]])[1:30,1:30] - t(temp3[1:30,1:30]) )^2 %>% as.vector)
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.NSBSR[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(temp3[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
})

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(PP.mat[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(PP.mat[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int.mat),4)
time.pp

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_PP(toy).Rdata"))

