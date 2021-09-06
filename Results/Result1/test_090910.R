

library(lattice)
library(dplyr)
rm(list=ls())
# temp = matrix(0,18,18)
for(nd in 1:100){
  
  nd <- 1
  
  setwd("D:/NSBSR_FINAL/THESIS_Results(New)/Result3/ANISOresults")
  load(file=paste0("NSBSR_Matern05_Grid(64,64)_ND=100_nd=",nd,".Rdata"))

  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  
  # thetahat = log(new.mcmcobj$Parset$lambda)
  # mat.theta = matrix(thetahat,n,n)
  thetahat  = Esti.result$statistics[(3+1+1):(3+1+145)]
  mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
  
  nh = ceiling(n/2)
  re.mat.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
                       cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
  dim(re.mat.theta)
  setwd("D:/NSBSR_FINAL")
  # pdf(paste0(Sys.Date(),"NSBSR-estimated_","IsoMatern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  # png(paste0(Sys.Date(),"NSBSR-estimated_","IsoMatern20_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".png"))
  setwd("D:/NSBSR_FINAL/THESIS_Results(New)/Result1")
  pdf(paste0(Sys.Date(),"NSBSR-estimated_AnIsoMatern05_Grid(32,32)_ND=1.pdf"))
  
  
  wireframe(re.mat.theta,default.scales =
              list(distance = 2,arrows = FALSE),shade=TRUE,zlab="",xlab="",ylab="",xlim=c(0,18),ylim=c(0,18))
  
  
  
  wireframe(re.mat.theta[seq(1,32,by=2),seq(1,32,by=2)],default.scales =
              list(distance = 2,arrows = FALSE),shade=TRUE,zlab="",xlab="",ylab="",xlim=c(0,32),ylim=c(0,32))
  
  wireframe((re.mat.theta*150+842),default.scales =
              list(distance = 1,arrows = FALSE),shade=TRUE,zlab="",xlab="",ylab="",xlim=c(0,18),ylim=c(0,18),zlim=c(-10,-3))
  setwd("D:/NSBSR_FINAL/THESIS_Results(New)/Result1")
  pdf(paste0(Sys.Date(),"NSBSR-estimated_IsoMatern05_Grid(64,64)_ND=1.pdf"))
  wireframe((re.mat.theta*180+1250),default.scales =
              list(distance = 1,arrows = FALSE),shade=TRUE,zlab="",xlab="",ylab="",xlim=c(0,32),ylim=c(0,32),zlim=c(-10,-3))
  dev.off()
  summary((re.mat.theta*180+1250)%>%as.vector)
  # 
  
}




rm(list=ls())
temp = matrix(0,18,18)
for(nd in 1:10){
  setwd("D:/NSBSR_FINAL/THESIS_Results(New)/Result3/ANISOresults")
  load(file=paste0("NSBSR_Matern20_Grid(32,32)_ND=100_nd=",nd,".Rdata"))
  
  Esti.result = summary(mcmc.list(NSBSR.EstiMat))
  
  # thetahat = log(new.mcmcobj$Parset$lambda)
  # mat.theta = matrix(thetahat,n,n)
  thetahat  = Esti.result$statistics[(3+1+1):(3+1+145)]
  mat.theta = matrix(c(thetahat[1:(length(thetahat)-1)],rev(thetahat)),n,n)
  
  nh = ceiling(n/2)
  re.mat.theta = rbind(cbind(mat.theta[rev(1:nh),rev(1:nh)],mat.theta[rev(1:nh),1:nh]),
                       cbind(mat.theta[1:nh,rev(1:nh)],mat.theta[1:nh,1:nh]))
  
  temp = temp + re.mat.theta
  # # pdf(paste0(Sys.Date(),"NSBSR-estimated_","AnIsoMatern05_","Grid(",grid.size[1],",",grid.size[2],")_ND=",nd,".pdf"),width=7*1,height=7*1)
  # par(mfrow=c(1,1),mar=c(1,1,1,1))
  # fields::image.plot(re.mat.theta,axes=F,cex.main=3.2,legend.mar=5.0)
  # # dev.off()
  # 
  
}
wireframe(temp+51,shade=TRUE,zlab="")
