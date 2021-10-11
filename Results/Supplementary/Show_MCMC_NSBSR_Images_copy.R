########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : Show_MCMC_NSBSR_Images.R

########################################################################################################
########################################################################################################

n1=simdata$D[1]; n2=simdata$D[2]; alpha=1/simdata$d[1]
gp2=gp(c(n1,n2),matern.specdens,c(0.5,alpha))
dens2=matern.specdens(gp2$omega,c(0.5,alpha),d=2)
dens2mat = matrix(dens2,n1,n2)
logdens2mat = log(dens2mat / sum(dens2mat))
nh = ceiling(n/2)
re.logdens2mat = rbind(cbind(logdens2mat[rev(1:nh),rev(1:nh)],logdens2mat[rev(1:nh),1:nh]),
                       cbind(logdens2mat[1:nh,rev(1:nh)],logdens2mat[1:nh,1:nh]))

# Settings...
iso = TRUE
n2 = simdata$nobs; mbar = ceiling(n2/2); n = sqrt(n2); m = ceiling((n+1)/2)
bdd = simdata$map$bdd; obs = simdata$map$obs; D = simdata$D
# data00 = data.frame(simdata$map,simdata$X,simdata$E,simdata$Y) %>% filter(obs==4)
data00 = data.frame(simdata$map,simdata$X,simdata$Y) %>% filter(obs==4)
loc = data00[,c("x","y")]-floor(D/2)

# Isotropic case
if( iso == TRUE ){
  
  truespd  <- matern.specdens(omega=loc,param=c((simdata$true.cov.par[3]/D[1]),simdata$true.cov.par[1]))
  truespd <- truespd / sum(truespd)
  
}
# AnIsotropic case
if( iso == FALSE ){
  
  # Rstar = matrix(c(-1/sqrt(5),1/sqrt(5),0,1),2,2)
  # Rstar = matrix(c(-1/sqrt(5),0,1/sqrt(5),1),2,2)
  Rstar = matrix(c(1,0,0,4),2,2)
  an.loc  = t(Rstar%*%t(loc))
  truespd  <- matern.specdens(omega=an.loc,param=c((simdata$true.cov.par[3]/D[1]),simdata$true.cov.par[1]))
  truespd <- truespd / sum(truespd)
  
}

#######################################################################
par(mfrow=c(2,3))
X = simdata$X; Y = simdata$Y

Yhat = PredResult$result1 + PredResult$result2
brk <- seq(min(Y,na.rm=T),max(Y,na.rm=T),length.out=64)
score = sqrt(mean((Y-Yhat)^2))

Ymat  = matrix(Y,   ncol=sqrt(simdata$nloc))
Yhmat = matrix(Yhat,ncol=sqrt(simdata$nloc))

col <- rev(rainbow(77)[1:63])

# Take 1
fields::image.plot(Ymat,main="TRUE")
# Take 2
fields::image.plot(Yhmat,main=paste0("MSE=",round(score,2)))
# Take 3
plot(Y,Yhat);abline(b=1,lty=2,col=2)
# Take 4
fields::image.plot(matrix(log(truespd),ncol=sqrt(length(truespd))),axes=F)

# Take 5
meantheta = Esti.summary$statistics[(new.mcmcobj$n.beta+1+1):(new.mcmcobj$n.beta+1+new.mcmcobj$n.theta),1]
if(n%%2==0){
  mat.theta  = matrix(c(meantheta,rev(meantheta)),ncol=n)
}else{
  mat.theta  = matrix(c(meantheta,rev(meantheta[1:(length(meantheta)-1)])),ncol=n)
} 
nh = ceiling(n/2)
rearranged.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
                         cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
rescaled.theta = -(-rearranged.theta-max(log(truespd)))^{1}*min(log(truespd))/min(-(-rearranged.theta-max(log(truespd)))^{1})
fields::image.plot(rescaled.theta,main="SPD(averaged)",axes=F)


# Take 6
lambda = new.mcmcobj$Parset$lambda
mat.lambda = matrix(lambda,ncol=n)/sum(lambda)
fields::image.plot(mat.lambda,main="exp(SPD)",axes=F)


#######################################################################
rm(bdd); rm(obs); rm(D)


