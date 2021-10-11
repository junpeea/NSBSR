
rm(list=ls())
setwd("~/Desktop/Code3(Acting)")
# load("200601StudyDataset.RData")
# load("2018AURAdataset.RData")
load("2017AURAdataset.RData")

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
# SPDE (INLA)
#######################################################################################################################
library(INLA) ## See http://www.r-inla.org/download/ for installation instructions.
library(fields) ## Only used for plotting
library(pals) ## Palettes for plotting

time.spde <- system.time({
  
  train.sim_data <- mysimdata[!is.na(mysimdata$O3L3),]; dim(train.sim_data)
  fit0 = lm(O3L3~R+S+NO2+VOC+UV+SO2,data=train.sim_data)
  train.sim_data$z = fit0$residuals
  train.sim_data$std <- 0.1
  head(train.sim_data); summary(train.sim_data)
  
  ## Model settings:
  spde.alpha <- 1.5 ## 1.5 gives an approximation of an exponential covariance
  
  ## inla() options:
  num.threads = 2 ## Limit memory usage by limiting the number of openmp threads.
  openmp.strategy = "large" ## "large" turns off nested parallelism
  
  n.row <- (xrange+1)
  n.col <- (yrange+1)
  indata <- mydata
  indata$z = NA
  indata$z[!is.na(mydata$O3L3)] <- train.sim_data$z
  extend <- -0.5
  max.edge <- 1/100*10
  
  ## Construct prior parameter distributions heuristically:
  ## Medium range, relative half-spread factor
  prior.range <- c(max(c(diff(range(indata$x)),
                         diff(unique(indata$y)))) * 1/5, 5)
  ## Medium sd, relative half-spread factor
  prior.sigma <- c(sd(indata$z, na.rm=TRUE) / 2, 4)
  
  ## Construct centred covariate versions
  LonCentre <- mean(range(indata$x))
  LatCentre <- mean(range(indata$y))
  indata$LonC <- indata$x - LonCentre
  indata$LatC <- indata$y - LatCentre
  
  lon <- matrix(indata$x,nrow=n.row)
  lat <- matrix(indata$y,nrow=n.row)
  temps <- list(mean=matrix(indata$O3L3,nrow=n.row))
  temps00 <- list(mean=matrix(indata$z,nrow=n.row))
  
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
      
      ok <- !is.na(indata$O3L3) | TRUE ## Include all; NA data are predicted
      data2 <- list(Temp=indata$z[ok], LonC=indata$LonC[ok], LatC=indata$LatC[ok],
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
                                    Temp=indata$O3L3*NA,
                                    SD=indata$O3L3*NA)
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
  setwd("~/Desktop/heatoncomparison/Code/SPDE")
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
  
  pred.SPDE  = matrix(X.pred%*%fit0$coefficients+pred2$mean,(xrange+1),(yrange+1))
  
})

Realdata = matrix(O3L3.full,(xrange+1),(yrange+1))
sdhatmat = matrix(pred2$sd,(xrange+1),(yrange+1))
mean((Realdata > pred.SPDE - 1.96 * sdhatmat) * (Realdata < pred.SPDE + 1.96 * sdhatmat))


par(mfrow=c(1,2))
fields::image.plot(O3L3.full,zlim=c(270,350))
fields::image.plot(pred.SPDE,zlim=c(270,350))

mae  = mean(abs((as.vector(O3L3.full) - as.vector(pred.SPDE))),na.rm=T);round(mae,3)
# mae  = mean(abs((as.vector(Krigemat) - as.vector(pred.SPDE))),na.rm=T);round(mae,3)
rmse = sqrt(mean((as.vector(O3L3.full) - as.vector(pred.SPDE))^2,na.rm=T));round(rmse,3)
# rmse = sqrt(mean((as.vector(Krigemat) - as.vector(pred.SPDE))^2,na.rm=T));round(rmse,3)
round(mean(int),3)
time.spde

setwd("~/Desktop/Code3(Acting)")
# save.image(paste0(Sys.Date(),"_SPDE.Rdata"))
# save.image(paste0(Sys.Date(),"_SPDE2018.Rdata"))
save.image(paste0(Sys.Date(),"_SPDE2017.Rdata"))

