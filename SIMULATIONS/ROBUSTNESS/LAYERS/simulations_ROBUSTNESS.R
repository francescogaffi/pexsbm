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
library(fastDummies)
library(greed)
library(randnet)
library(Rcpp)
library(RcppArmadillo)
library(SpecClustPack)

# remove all except for functions
rm(list = setdiff(ls(), lsf.str()))

# load source functions
source("pex-sbm.R")
source("casc.R")

##########################################################################
################# SIMULATION ROBUSTNESS UNINFORMATIVE LAYERS #############
##########################################################################
load("simulations_perm.Rdata")

# setting the dimensions and other useful quantities
n_iter  <- 10000
burn_in <- 2000
V <- 80
N  <- c(30,30,15,5)

# for scenario 3 (uninformative layers)

# for CASC
X_casc <- as.matrix(dummy_cols(lay_perm)[,-1])

# HPD hyperparameters setting
theta  <- 0.5
theta0 <- 4
sigma  <- 0.2
sigma0 <- 0.8
alpha   <- 5
lambda  <- 10
alpha0  <- 12
lambda0 <- 3

# define the matrices containing the performance measures
vi            <- matrix(0,4,10)
z_estimated   <- array(0,c(4,10,80))


for (r in 1:10){
	
Y1 <- Y1_replicates[,,r]
Y3 <- Y1[sel,sel]

#---------------------------------------------------------------
# fit pEx-SBM(HDP) to scenario 3
set.seed(r)
model3       <- pex_sbm(Y3,n_iter,N,1,1,theta,theta0)
z_post3      <- model3[[1]]
P_post3      <- pr_cc(z_post3)
z_hat3       <- minVI(P_post3,method="avg",max.k=20)$cl

vi[1,r] <- VI(z_perm,t(z_hat3))
z_estimated[1,r,] <- z_hat3

#---------------------------------------------------------------
# fit pEx-SBM(HNSP) to scenario 3
set.seed(r)
model3       <- pex_sbm(Y3,n_iter,N,1,1,sigma=sigma, sigma0=sigma0, prior="HNSP")
z_post3      <- model3[[1]]
P_post3      <- pr_cc(z_post3)
z_hat3       <- minVI(P_post3,method="avg",max.k=20)$cl

vi[2,r] <- VI(z_perm,t(z_hat3))
z_estimated[2,r,] <- z_hat3

#---------------------------------------------------------------
# fit pEx-SBM(HDP hyperprior) to scenario 3
set.seed(r)
model3       <- hdpgamma_sbm(Y3,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post3      <- model3[[1]]
P_post3      <- pr_cc(z_post3)
z_hat3       <- minVI(P_post3,method="avg",max.k=20)$cl

vi[3,r] <- VI(z_perm,t(z_hat3))
z_estimated[3,r,] <- z_hat3

#---------------------------------------------------------------
# CASC scenario 3 - Binkiewicz et al. (2017)
K <- 8
set.seed(r)
cluster_casc <- casc(Y3,X_casc,K)$cluster

vi[4,r] <- VI(z_perm,t(cluster_casc))
z_estimated[4,r,] <- cluster_casc

print(r)}

save(vi, z_estimated, file="output_ROBUSTNESS.RData")
