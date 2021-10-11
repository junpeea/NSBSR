########################################################################################################
########################################################################################################
# Thesis : Bayesian spatial prediction with nonparametric modeling of a spectral density

# PROJECT name : NSBSR

# Source name : Get2dGridEnvironment.R

########################################################################################################
########################################################################################################
# try(if(sum(grid.size/obs.delta==floor(grid.size/obs.delta))!=2) stop("ERROR: Invalid resolutions settings"))

# sp1 = seq(0,grid.size[1],by=1); sp2 = seq(0,grid.size[2],by=1)
# sp1 = seq(1,grid.size[1],by=1); sp2 = seq(1,grid.size[2],by=1)
sp1 = seq(0,(grid.size[1]-1),by=1); sp2 = seq(0,(grid.size[2]-1),by=1)
space.grid = expand.grid(sp1,sp2)
omega.grid = expand.grid(2*pi*sp1/grid.size[1],2*pi*sp2/grid.size[2])

nloc = nrow(space.grid)

bddind = rep(1,nloc)
bddind[which(space.grid[,1] %% (grid.size[1]/2) == 0 )] <- 2
bddind[which(space.grid[,2] %% (grid.size[2]/2) == 0 )] <- 2
plot(space.grid,col=bddind,main="Show_Boundaries")

obsind = rep(NA,nloc)
obsind[which(rowSums(space.grid %% obs.delta) == 0 )] <- 4
nobs = length(which(is.na(obsind)==0))
set.seed(1); p = rbinom(nobs,1,obs.prob)
obsind[which(obsind==4)] <- 4*p + 2*(1-p)
plot(space.grid,col=obsind,main="Show_MonitoredLocation")

map.grid = data.frame(space.grid,omega.grid,bddind,obsind); colnames(map.grid) <- c("x","y","w.x","w.y","bdd","obs")

d.o.o.w   <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.x","w.y")) %>% dist %>% as.matrix
d.o.o.w.x <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.x")) %>% dist %>% as.matrix
d.o.o.w.y <- map.grid %>% filter(obs > 0) %>% dplyr::select(c("w.y")) %>% dist %>% as.matrix



########################################################################################################
rm(sp1);rm(sp2);rm(obsind);rm(bddind)
