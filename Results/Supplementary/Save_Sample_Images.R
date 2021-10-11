########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : Save_Sample_Images.R

########################################################################################################
########################################################################################################

filename = paste0("Dataset/Isotropic/Nugget_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","Nugget",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
        
  }
  cat(date(),"Isotropic Nugget Sample Images updated","\n")
}

filename = paste0("Dataset/Isotropic/Spheric_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","Spheric",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
    
  }
  cat(date(),"Isotropic Shperic Sample Images updated","\n")
}

filename = paste0("Dataset/Isotropic/Matern01_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","Matern01",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
    
  }
  cat(date(),"Isotropic Matern(Bumpy) Sample Images updated","\n")
}

filename = paste0("Dataset/Isotropic/Matern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","Matern05",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
    
  }
  cat(date(),"Isotropic Exponential Sample Images updated","\n")
}

filename = paste0("Dataset/Isotropic/Matern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","Matern20",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
    
  }
  cat(date(),"Isotropic Matern(Smooth) Sample Images updated","\n")
}

filename = paste0("Dataset/Isotropic/SqExp05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","SqExp05",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
    
  }
  cat(date(),"Isotropic Squared Exponential(Bumpy) Sample Images updated","\n")
}

filename = paste0("Dataset/Isotropic/SqExp15_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","SqExp15",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
    
  }
  cat(date(),"Isotropic Squared Exponential(Smooth) Sample Images updated","\n")
}

filename = paste0("Dataset/Isotropic/Gauss_","Grid(",grid.size[1],",",grid.size[2],")_ND=",ND,".Rdata")
if( file.exists(filename) ){
  
  load(file=filename)
  
  for(nd in 1:ND){
    
    pdf(paste("Suppl/Sample_Images/Isotropic","Gauss",nd,"image.pdf",sep="_"))
    
    fields::image.plot(matrix(Y[[nd]],nrow=sqrt(nloc)),main=paste0("Show","_Y_","ND=",nd))
    
    dev.off()
    
  }
  cat(date(),"Isotropic Gauss Sample Images updated","\n")
}


########################################################################################################
rm(filename)