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

load("simulations.Rdata")

# setting the dimensions and other useful quantities
n_iter  <- 10000
burn_in <- 2000
V <- 80
N  <- c(30,30,15,5)

# HPD hyperprior parameters setting
alpha   <- 5
lambda  <- 10
alpha0  <- 12
lambda0 <- 3

# define the matrices containing the performance measures
waic          <- matrix(0,2,10)
vi            <- matrix(0,2,10)
expected_vi   <- matrix(0,2,10)
vi_bound      <- matrix(0,2,10)
z_estimated   <- array(0,c(2,10,80))
H_estimated   <- matrix(0,2,10)
H_posterior   <- array(0,c(2,10,3))
misclass      <- matrix(0,1,10)
pred_prob_err <- matrix(0,1,10)

for (r in 1:10){
	
Y1 <- Y1_replicates[,,r]
Y2 <- Y2_replicates[,,r]

#---------------------------------------------------------------
# fit pEx-SBM(HDP hyperprior) to scenario 1
set.seed(r)
model1       <- hdpgamma_sbm(Y1,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post1      <- model1[[1]]
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
# fit pEx-SBM(HDP hyperprior) to scenario 2
set.seed(r)
model2       <- hdpgamma_sbm(Y2,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post2      <- model2[[1]]
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

#---------------------------------------------------------------
# group prediction (scenario 1)
# remove m nodes at random
set.seed(r)
m      <- 10
Vrem   <- V-m
rem    <- randsample(V,m)
Yrem   <- Y1[-rem,-rem]
layrem <- lay[-rem]
remlay <- lay[rem]
Nrem   <- as.numeric(table(layrem))

# fit pEx-SBM(HDP hyperprior)
set.seed(r)
model_rem <- hdpgamma_sbm(Yrem,n_iter,Nrem,1,1,alpha,lambda,alpha0,lambda0)
z_post    <- model_rem[[1]]
w_post    <- model_rem[[2]]
theta_post  <- model_rem[[3]]
theta0_post <- model_rem[[4]]

# initialize co-clustering matrix frames
  P_new     <- matrix(0,length(rem),length(rem))
  P_new_sum <- matrix(0,length(rem),length(rem))
  P_old     <- matrix(0,length(rem),sum(Nrem))
  P_old_sum <- matrix(0,length(rem),sum(Nrem))

# mc iterations
n_mc <- n_iter - burn_in

# for each posterior allocation sample
for (i in 1:n_mc){
    n <- i+burn_in
    #get matrices l and q
    summary <- getsummary(z_post[,n],w_post[,n],layrem)
    l       <- summary[[1]]
    q       <- summary[[2]]
  
  for(ind1 in 1:m){
    # compute co-clustering probabilities with old nodes
      for(ind2 in 1:Vrem){
        P_old[ind1,ind2] <- coclust_old(l,q,remlay[ind1],z_post[ind2,n],theta_post[n],theta0_post[n],NA,NA,"HDP",Nrem)
      }
    # compute co-clustering probabilities among new nodes
      if(ind1>1)
        for(ind3 in 1:(ind1-1)){
          P_new[ind1,ind3] <- coclust_new(l,q,remlay[ind1],remlay[ind3],theta_post[n],theta0_post[n],NA,NA,"HDP",Nrem)
        }
  }
  
    P_new_sum <- P_new_sum+P_new
    P_old_sum <- P_old_sum+P_old
}
  P_new       <- P_new_sum/n_mc
  P_new       <- P_new+t(P_new)
  diag(P_new) <- rep(1,m)
  P_old       <- P_old_sum/n_mc

# frame the co-clustering matrix of old nodes with co-clustering probabilities
  P <- rbind(pr_cc(z_post[,(burn_in+1):n_iter]),P_old)
  P <- cbind(P,rbind(t(P_old),P_new))

z_hat_prediction  <- minVI(P,method="avg",max.k=20)$cl
misclass[1,r] <- sum(z[rem]!=z_hat_prediction[71:80])

#---------------------------------------------------------------
# edge prediction
  set.seed(r)
  pred      <- pred_conn_hyp(Yrem,layrem,remlay,model_rem,n_mc)
  Ynew      <- pred[[1]]
  perm      <- c(c(1:80)[-rem],rem)

# compare with adjacency matrix of old nodes framed with true connection probabilities
Phi_true <- Phi[z[perm],z[perm]]
pred_prob_err[1,r]<-(sum((Ynew[71:80,1:70]-Phi_true[71:80,1:70])^2)+sum(lowerTriangle((Ynew[71:80,71:80]-Phi_true[71:80,71:80])^2)))/(length(c(Ynew[71:80,1:70]))+length(c(lowerTriangle(Ynew[71:80,71:80]))))

print(r)}

save(waic, vi, expected_vi, vi_bound, z_estimated, H_estimated, H_posterior, misclass, pred_prob_err, file="output_HDP_hyperprior.RData")
