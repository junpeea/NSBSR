#######################################################################################################################
# SPDE (INLA)
#######################################################################################################################

time.spde <- system.time({
  n <- 30
  temp = expand.grid(c(0:(n-1)),c(0:(n-1)))
  sim_process = data.frame(x=temp$Var1,y=temp$Var2)
  sim_process$Y = as.vector(Y0)
  
  ## Model settings:
  spde.alpha <- 1.5 ## 1.5 gives an approximation of an exponential covariance
  
  ## inla() options:
  num.threads = 2 ## Limit memory usage by limiting the number of openmp threads.
  openmp.strategy = "large" ## "large" turns off nested parallelism
  
  n.row <- 30
  n.col <- 30
  indata <-sim_process
  indata$Y00  <- as.vector(Y00)
  extend <- -0.5
  max.edge <- 1/100*10
  
  ## Construct prior parameter distributions heuristically:
  ## Medium range, relative half-spread factor
  prior.range <- c(max(c(diff(range(indata$x)),
                         diff(unique(indata$y)))) * 1/5, 5)
  ## Medium sd, relative half-spread factor
  prior.sigma <- c(sd(indata$Y00, na.rm=TRUE) / 2, 4)
  
  ## Construct centred covariate versions
  LonCentre <- mean(range(indata$x))
  LatCentre <- mean(range(indata$y))
  indata$LonC <- indata$x - LonCentre
  indata$LatC <- indata$y - LatCentre
  
  lon <- matrix(indata$x,nrow=n.row)
  lat <- matrix(indata$y,nrow=n.row)
  temps <- list(mean=matrix(indata$Y,nrow=n.row))
  temps00 <- list(mean=matrix(indata$Y00,nrow=n.row))
  
  # install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
  library(INLA) ## See http://www.r-inla.org/download/ for installation instructions.
  library(fields) ## Only used for plotting
  library(pals) ## Palettes for plotting
  
  
  ## Run the full model.
  time2 <- system.time({
    time2.prep <- system.time({                                                                                            
      loc2 <- cbind(indata$x, indata$y)
      mesh2 <- inla.mesh.create(loc=loc2,
                                extend=list(offset=extend),
                                refine=list(min.angle=30))
      proj2 <- inla.mesh.projector(mesh2, loc2)
      col.idx <- apply(proj2$proj$bary, 1, which.max)
      idx2 <- mesh2$graph$tv[proj2$proj$t + (col.idx-1)*nrow(mesh2$graph$tv)]
      
      ## lognormal prior
      ## Centre the parameterisation at range=1, sigma=1
      nu <- spde.alpha - 1
      kappa.zero <- sqrt(8*nu) / 1
      tau.zero <- (gamma(nu) / (gamma(spde.alpha) * 4*pi * kappa.zero^(2*nu)) )^0.5 / 1
      ## sigma^2 = Gamma(0.5)/( Gamma(1.5) (4\pi)^(dim/2) * kappa^(2*0.5) * tau^2 )
      ## tau = [ Gamma(0.5)/( Gamma(1.5) (4\pi)^(dim/2) * kappa^(2*0.5) ) ]^0.5 / sigma
      ## kappa = sqrt(8*0.5) / range
      ## tau = [ Gamma(0.5)/( Gamma(1.5) (4\pi)^(dim/2) * (sqrt(8*0.5)/range)^(2*0.5) ) ]^0.5 / sigma
      B.tau   <- cbind(log(tau.zero), nu, -1)
      B.kappa <- cbind(log(kappa.zero), -1, 0)
      
      theta.prior.mean <- log(c(prior.range[1], prior.sigma[2]))
      ## Half-width of prior prediction interval, on log-scale: 2*sd = log(rel)
      theta.prior.prec <- 4 / log(c(prior.range[2], prior.sigma[2]))^2
      spde2 <- inla.spde2.matern(mesh2, alpha=spde.alpha,
                                 B.tau = B.tau, B.kappa = B.kappa,
                                 theta.prior.mean = theta.prior.mean,
                                 theta.prior.prec = theta.prior.prec,
                                 constr=FALSE)
      
      hyper.range.initial <- log(prior.range[1])
      hyper.sigma.initial <- log(prior.sigma[1])
      hyper.family.initial <- 2
      
      spde2$f$hyper.default$theta1$initial <- hyper.range.initial
      spde2$f$hyper.default$theta2$initial <- hyper.sigma.initial
      
      formula2 <- Temp ~ 1 + LonC + LatC + f(field, model=spde)
      
      ok <- !is.na(indata$Y) | TRUE ## Include all; NA data are predicted
      data2 <- list(Temp=indata$Y[ok], LonC=indata$LonC[ok], LatC=indata$LatC[ok],
                    field=idx2[ok], spde=spde2)
    })
    time2.run <- system.time({  
      result2 <- inla(formula2, data=data2, family="gaussian",
                      control.family=list(hyper=list(theta=list(initial=hyper.family.initial))),
                      control.predictor=list(compute=TRUE),
                      control.inla=list(strategy="gaussian", int.strategy="eb",
                                        force.diagonal=TRUE, stupid.search=FALSE),
                      verbose=TRUE, num.threads=num.threads,
                      control.compute=list(openmp.strategy=openmp.strategy))
    })
    times2.post <- system.time({
      reconstruction2 <- data.frame(Lon=indata$x,
                                    Lat=indata$y,
                                    Temp=indata$Y*NA,
                                    SD=indata$Y*NA)
      reconstruction2$Temp[ok] <- result2$summary.linear.predictor[,"mean"]
      reconstruction2$SD[ok] <- 
        sqrt(result2$summary.linear.predictor[,"sd"]^2
             + 1/result2$summary.hyperpar[1,"0.5quant"])
      
      temps2 <- list(mean=matrix(reconstruction2$Temp,nrow=n.row),
                     sd=matrix(reconstruction2$SD,nrow=n.row))
      
      ## Convert model parameters to humanly readable form
      result2.field <- inla.spde.result(result2, "field", spde=spde2)
      param2 <- rbind("nugget.sd"=vapply(2:7, function(x) {
        if (x==1) {
          1
        } else if (x==8) {
          0
        } else if (x==2) {
          inla.emarginal(function(x) x^(-0.5),
                         result2$marginals.hyperpar[[1]])
        } else if (x==3) {
          m <- inla.emarginal(function(x) x^(-0.5),
                              result2$marginals.hyperpar[[1]])
          inla.emarginal(function(x) (x^(-0.5) - m)^2,
                         result2$marginals.hyperpar[[1]])^0.5
        } else {
          kk <- c(5,4,3,6)[x-3]
          result2$summary.hyperpar[1,kk]^(-0.5)
        }
      }, 1.0))
      param2 <- rbind(param2, "field.range"=vapply(2:7, function(k) {
        if (k==1 | k==8) {
          result2.field$summary.log.range.nominal[1, k]
        } else if (k==2) {
          inla.emarginal(function(x) x,
                         result2.field$marginals.range.nominal[[1]])
        } else if (k==3) {
          m <- inla.emarginal(function(x) x,
                              result2.field$marginals.range.nominal[[1]])
          inla.emarginal(function(x) (x - m)^2,
                         result2.field$marginals.range.nominal[[1]])^0.5
        } else {
          exp(result2.field$summary.log.range.nominal[,k])
        }
      }, 1.0))
      param2 <- rbind(param2,
                      "field.sd"=vapply(2:7, function(k) {
                        if (k==1 | k==8) {
                          result2.field$summary.log.variance.nominal[1, k]
                        } else if (k==2) {
                          inla.emarginal(function(x) x^0.5, result2.field$marginals.variance.nominal[[1]])
                        } else if (k==3) {
                          m <- inla.emarginal(function(x) x^0.5,
                                              result2.field$marginals.variance.nominal[[1]])
                          inla.emarginal(function(x) (x^0.5 - m)^2,
                                         result2.field$marginals.variance.nominal[[1]])^0.5
                        } else {
                          exp(result2.field$summary.log.variance.nominal[,k] / 2)
                        }
                      }, 1.0))
      colnames(param2) <- colnames(result2.field$summary.log.range.nominal)[-c(1,8)]
      param2 <- as.data.frame(param2)
    })
  })
  
  
  #### Assess predictions #####
  setwd("C:/Users/yjun/Documents/GitHub/heatoncomparison/Code/SPDE")
  source("assessment.R")
  
  calc.scores <- function(pred, y, name) {
    s <- data.frame(c(
      MSE=mean(sqerr(pred, y)),
      MAE=mean(abserr(pred, y)),
      IGN=mean(ign(pred, y)),
      CRPS=mean(crps(pred, y)),
      Interval=mean(inter(pred, y)),
      Coverage=mean(cover(pred, y))))
    colnames(s) <- name
    s
  }
  
  pred2 <- pred.obj(reconstruction2$Temp, reconstruction2$SD)
  int <- 2*1.96*pred2$sd 
  
  sdhatmat = matrix(pred2$sd,30,30)
  mean((Realdata > pred2$mean - 1.96 * sdhatmat) * (Realdata < pred2$mean + 1.96 * sdhatmat))
  
  
  mx  = mean(sim_process$x[which(ind==FALSE)]) ;my  = mean(sim_process$y[which(ind==FALSE)])
  sdx = sd(sim_process$x[which(ind==FALSE)])   ;sdy = sd(sim_process$y[which(ind==FALSE)])
  X = as.matrix(cbind(1,(sim_process$x-mx)/sdx,(sim_process$y-my)/sdy))
  ind = is.na(as.vector(Y00));summary(ind)
  sim_data <- sim_process[which(ind==FALSE),]
  sim_data$x1 = rep(1,nrow(sim_data))
  sim_data$x2 = scale(sim_process$x[which(ind==FALSE)])
  sim_data$x3 = scale(sim_process$y[which(ind==FALSE)])
  fit0 = lm(Y~-1+x1+x2+x3,data=sim_data)
  pred.SPDE  = matrix(pred2$mean,30,30)
  pred.NSBSR = matrix(X%*%fit0$coefficients + as.vector(temp2),30,30)
  
  summary((t(mylist[[1]])[1:30,1:30] - t(pred.SPDE[1:30,1:30]) )^2 %>% as.vector)
  cor(x=t(mylist[[1]])[1:30,1:30]%>%as.vector,y=t(pred.SPDE[1:30,1:30])%>%as.vector,use="complete.obs")
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.SPDE[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
  summary((t(mylist[[1]])[1:30,1:30] - t(pred.NSBSR[1:30,1:30]))^2 %>% as.vector)
  cor(x=t(mylist[[1]])[1:30,1:30]%>%as.vector,y=t(pred.NSBSR[1:30,1:30])%>%as.vector,use="complete.obs")
  # summary((t(mylist[[1]])[1:30,1:30] - t(temp3[1:30,1:30]) )^2 %>% as.vector)
  par(mfrow=c(2,2))
  fields::image.plot(Y0[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y00[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(Y0[seq(1,30,by=2),seq(1,30,by=2)],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  fields::image.plot(pred.NSBSR[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  # fields::image.plot(temp3[1:29,1:29],zlim=c(0.99*min(Y0,na.rm=T),1.01*max(Y0,na.rm=T)))
  
})

mae  = mean(abs((t(mylist[[k]])[1:30,1:30] - t(pred.SPDE[1:30,1:30]))),na.rm=T);round(mae,4)
rmse = sqrt(mean((t(mylist[[k]])[1:30,1:30] - t(pred.SPDE[1:30,1:30]) )^2,na.rm=T));round(rmse,4)
round(mean(int),4)
time.spde

setwd("C:/Users/yjun/Desktop/WORK2021/NSBSR_FINAL_210702/MODIS")
save.image(paste0(Sys.Date(),"_SPDE(toy).Rdata"))



