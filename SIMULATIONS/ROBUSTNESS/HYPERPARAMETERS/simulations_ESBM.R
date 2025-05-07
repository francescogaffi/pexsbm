# load useful libraries
library(RColorBrewer)
library(igraph)
library(pheatmap)
library(gridExtra)
library(grid)
library(copula)
library(mcclust.ext)
library(reshape)
library(gdata)
library(igraph)
library(lattice)
library(cowplot)
library(ggplot2)
library(coda)
library(greed)
library(LaplacesDemon)
library(tidyverse)
library(pracma)
library(dummy)

# remove all except for functions
rm(list = setdiff(ls(), lsf.str()))

# load source functions
source("pex-sbm.R")
source("numcl.R")

##########################################################################
################# SIMULATION SCENARIOS FOR ESBM DP #######################
##########################################################################

load("simulations.Rdata")

# setting the dimensions and other useful quantities
n_iter  <- 10000
burn_in <- 2000
V <- 80
N  <- c(30,30,15,5)

# ESBM hyperparameters setting
theta_esbm <- 2.75
expnumcl_dp(V,theta_esbm)

# define the matrices containing the performance measures
waic          <- matrix(0,2,10)
vi            <- matrix(0,2,10)
expected_vi   <- matrix(0,2,10)
vi_bound      <- matrix(0,2,10)
z_estimated   <- array(0,c(2,10,80))
H_estimated   <- matrix(0,2,10)
H_posterior   <- array(0,c(2,10,3))

for (r in 1:10){
	
Y1 <- Y1_replicates[,,r]
Y2 <- Y2_replicates[,,r]

#---------------------------------------------------------------
# fit ESBM to scenario 1
model1       <- esbm(Y1,r,n_iter,"DP",alpha_PY = theta_esbm, sigma_PY = 0,x=lay, alpha_xi = rep(1,4))
z_post1      <- model1
P_post1      <- pr_cc(z_post1)
z_hat1       <- minVI(P_post1,method="avg",max.k=20)$cl

vi_bound[1,r] <- credibleball(z_hat1,t(z_post1[,(burn_in+1):n_iter]))[[5]]
expected_vi[1,r] <- VI(z,t(z_post1[,(burn_in+1):n_iter]))
vi[1,r] <- VI(z,t(z_hat1))
H_estimated[1,r] <- max(z_hat1)
H_posterior[1,r,] <- summary(apply(z_post1[,(burn_in+1):n_iter],2,max))[c(2,3,5)]
z_estimated[1,r,] <- z_hat1

LL <- matrix(nrow=V*(V-1)/2,ncol=n_iter-burn_in)
z_waic <- z_post1[,(burn_in+1):n_iter]
for (t in 1:dim(z_waic)[2]){
  LL[,t]<-sampleLL(z_waic[,t],Y1,1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
waic[1,r] <- WAIC(LL)$WAIC

#---------------------------------------------------------------
# fit ESBM to scenario 2
set.seed(r)
model2       <- esbm(Y2,r,n_iter,"DP",alpha_PY = theta_esbm, sigma_PY = 0,x=lay, alpha_xi = rep(1,4))
z_post2      <- model2
P_post2      <- pr_cc(z_post2)
z_hat2       <- minVI(P_post2,method="avg",max.k=20)$cl

vi_bound[2,r] <- credibleball(z_hat2,t(z_post2[,(burn_in+1):n_iter]))[[5]]
expected_vi[2,r] <- VI(z,t(z_post2[,(burn_in+1):n_iter]))
vi[2,r] <- VI(z,t(z_hat2))
H_estimated[2,r] <- max(z_hat2)
H_posterior[2,r,] <- summary(apply(z_post2[,(burn_in+1):n_iter],2,max))[c(2,3,5)]
z_estimated[2,r,] <- z_hat2

LL <- matrix(nrow=V*(V-1)/2,ncol=n_iter-burn_in)
z_waic <- z_post2[,(burn_in+1):n_iter]
for (t in 1:dim(z_waic)[2]){
  LL[,t]<-sampleLL(z_waic[,t],Y2,1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
waic[2,r] <- WAIC(LL)$WAIC

print(r)}

save(waic, vi, expected_vi, vi_bound, z_estimated, H_estimated, H_posterior, file="output_ESBM.RData")
