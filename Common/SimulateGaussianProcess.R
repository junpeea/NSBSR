########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : SimulateGaussianProcess.R

########################################################################################################
########################################################################################################

seedvec = Y.seed*c(1:ND)
beta = c(0.01,0.02,0.03); sig2_eps = 1; cor.par=c(0.5,(grid.size[1]+grid.size[2])/2.5)


filename = paste0("Dataset/Isotropic/Nugget_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  # model.option = RMnugget(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  model.option = RMnugget(var=sig2_eps,Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"Isotropic Nugget Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/Nugget_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=0.25)
  # model.option = RMnugget(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  model.option = RMnugget(var=sig2_eps,Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"AnIsotropic Nugget Sample Dataset generated","\n")
}


filename = paste0("Dataset/Isotropic/Spheric_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  model.option = RMspheric(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"Isotropic Spheric Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/Spheric_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=0.25)
  model.option = RMspheric(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"AnIsotropic Spheric Sample Dataset generated","\n")
}


filename = paste0("Dataset/Isotropic/Matern01_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  cor.par[1] = 0.1
  model.option = RMmatern(nu=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"Isotropic Matern(Bumpy) Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/Matern01_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=0.25)
  cor.par[1] = 0.1
  model.option = RMmatern(nu=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"AnIsotropic Matern(Bumpy) Sample Dataset generated","\n")
}


filename = paste0("Dataset/Isotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  cor.par[1] = 0.5
  model.option = RMmatern(nu=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"Isotropic Exponential Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=4)
  cor.par[1] = 0.5
  model.option = RMmatern(nu=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"AnIsotropic Exponential Sample Dataset generated","\n")
}


filename = paste0("Dataset/Isotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  cor.par[1] = 2.0
  model.option = RMmatern(nu=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option); save.image(file=filename)
  
  cat(date(),"Isotropic Matern(Smooth) Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=pi/4,ratio=5)
  cor.par[1] = 2.0
  model.option = RMmatern(nu=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option); save.image(file=filename)
  
  cat(date(),"AnIsotropic Matern(Smooth) Sample Dataset generated","\n")
}


filename = paste0("Dataset/Isotropic/SqExp05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  cor.par[1] = 0.5
  model.option = RMpoweredexp(alpha=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"Isotropic Squared Exponential(Bumpy) Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/SqExp05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=0.25)
  cor.par[1] = 0.5
  model.option = RMpoweredexp(alpha=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"AnIsotropic Squared Exponential(Bumpy) Sample Dataset generated","\n")
}


filename = paste0("Dataset/Isotropic/SqExp15_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  cor.par[1] = 1.5
  model.option = RMpoweredexp(alpha=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"Isotropic Squared Exponential(Smooth) Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/SqExp15_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=0.25)
  cor.par[1] = 1.5
  model.option = RMpoweredexp(alpha=cor.par[1],var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"AnIsotropic Squared Exponential(Smooth) Sample Dataset generated","\n")
}


filename = paste0("Dataset/Isotropic/Gauss_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=1)
  model.option = RMgauss(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"Isotropic Gauss Sample Dataset generated","\n")
}
filename = paste0("Dataset/AnIsotropic/Gauss_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( !file.exists(filename) ){
  
  aniso.option = RMangle(angle=0,ratio=0.25)
  model.option = RMgauss(var=sig2_eps,scale=cor.par[2],Aniso=aniso.option)
  
  El = seedvec %>% lapply(function(seed){
    set.seed(seed);epsi = RandomFields::RFsimulate(model=model.option,x=map.grid$x,y=map.grid$y)$variable1;return(epsi)
  })
  Yl = seedvec %>% lapply(function(w) (X %>% as.matrix)%*%beta + sqrt(sig2_eps) * El[[w]])
  rm(model.option);save.image(file=filename)
  
  cat(date(),"AnIsotropic Gauss Sample Dataset generated","\n")
}

########################################################################################################
rm(filename)



