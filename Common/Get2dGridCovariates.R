########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : Get2dGridCovariates.R

########################################################################################################
########################################################################################################

seed = X.seed

X0 = rep(1,nloc)
set.seed(seed); p = rbinom(nloc,1,0.5); seed=seed+123
set.seed(seed); X1 = p*rnorm(nloc,0,1) + (1-p)*rnorm(nloc,0,5); seed=seed+123
set.seed(seed); X2 = rexp(nloc,1)
X = data.frame(X0,X1,X2)

fields::image.plot(matrix(X0,nrow=sqrt(nloc)),main="Image_X0")
fields::image.plot(matrix(X1,nrow=sqrt(nloc)),main="Image_X1")
fields::image.plot(matrix(X2,nrow=sqrt(nloc)),main="Image_X2")

########################################################################################################
rm(seed);rm(X.seed);rm(p)