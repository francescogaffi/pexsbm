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
####### SIMULATION SCENARIOS FOR PEX-SBM WITH HPD HYPERPRIOR #############
##########################################################################

load("simulations_consistency.Rdata")

# setting the dimensions and other useful quantities
n_iter  <- 10000
burn_in <- 2000

# HPD hyperprior parameters setting
alpha   <- 5
lambda  <- 10
alpha0  <- 12
lambda0 <- 3

# define the matrices containing the performance measures
vi            <- array(0,c(2,10,7))

for (r in 1:7){
	
Y1_temp <- Y1_consistency[[r]]
Y2_temp <- Y2_consistency[[r]]

lay <- lay_consistency[[r]]
z <- z_consistency[[r]]
N <- c(table(lay))

for (s in 1:10){
	
Y1 <- Y1_temp[,,s]
Y2 <- Y2_temp[,,s]
V <- dim(Y1)[1]

#---------------------------------------------------------------
# fit pEx-SBM(HDP hyperprior) to scenario 1
set.seed(s)
model1       <- hdpgamma_sbm(Y1,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post1      <- model1[[1]]
P_post1      <- pr_cc(z_post1)
z_hat1       <- minVI(P_post1,method="avg",max.k=20)$cl

vi[1,s,r] <- VI(z,t(z_hat1))


#---------------------------------------------------------------
# fit pEx-SBM(HDP hyperprior) to scenario 2
set.seed(s)
model2       <- hdpgamma_sbm(Y2,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post2      <- model2[[1]]
P_post2      <- pr_cc(z_post2)
z_hat2       <- minVI(P_post2,method="avg",max.k=20)$cl

vi[2,s,r] <- VI(z,t(z_hat2))

print(s)}
print(r)}

save(vi, file="output_HDP_hyperprior.RData")