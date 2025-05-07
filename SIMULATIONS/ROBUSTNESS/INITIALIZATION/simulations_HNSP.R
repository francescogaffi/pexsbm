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

##########################################################################
################# SIMULATION SCENARIOS FOR PEX-SBM WITH HNSP #############
##########################################################################

load("simulations.Rdata")

# setting the dimensions and other useful quantities
n_iter  <- 10000
burn_in <- 2000
V <- 80
N  <- c(30,30,15,5)

# HNSP hyperparameters setting
sigma  <- 0.2
sigma0 <- 0.8

# define the matrices containing the performance measures
expected_vi   <- matrix(0,2,10)

for (r in 1:10){
	
Y1 <- Y1_replicates[,,r]
Y2 <- Y2_replicates[,,r]

#---------------------------------------------------------------
# fit pEx-SBM(HNSP) to scenario 1
set.seed(r)
model1       <- pex_sbm(Y1,n_iter,N,1,1,sigma=sigma, sigma0=sigma0, prior="HNSP",z_init = rep(1,V))
z_post1      <- model1[[1]]
P_post1      <- pr_cc(z_post1)
z_hat1       <- minVI(P_post1,method="avg",max.k=20)$cl

expected_vi[1,r] <- VI(z,t(z_post1[,(burn_in+1):n_iter]))

#---------------------------------------------------------------
# fit pEx-SBM(HNSP) to scenario 2
set.seed(r)
model2       <- pex_sbm(Y2,n_iter,N,1,1,sigma=sigma, sigma0=sigma0, prior="HNSP",z_init = rep(1,V))
z_post2      <- model2[[1]]
P_post2      <- pr_cc(z_post2)
z_hat2       <- minVI(P_post2,method="avg",max.k=20)$cl

expected_vi[2,r] <- VI(z,t(z_post2[,(burn_in+1):n_iter]))

print(r)}

save(expected_vi,  file="output_HNSP.RData")
