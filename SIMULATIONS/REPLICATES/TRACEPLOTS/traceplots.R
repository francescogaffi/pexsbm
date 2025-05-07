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

#-----------------------------------------------------------------------------------
# FUNCTION TO COMPUTE BETA-BINOMIAL LIKELIHOOD  ------------------------------------
#-----------------------------------------------------------------------------------

log_pY_z <- function(Y,z,a,b){
# in: Adjacency matrix Y, one vector of node labels z, hyperparameters (a,b) of Beta priors for block probabilities
# out: logarithm of the beta-binomial likelihood for (Y|z) 
options(warn=-1)

H <- length(unique(z))
colnames(Y) <- rownames(Y) <- z

edge_counts <- melt(Y)

Y_c <- 1 - Y
diag(Y_c) <- 0
non_edge_counts <- melt(Y_c)
	
Edge <- matrix(aggregate(edge_counts[,3],by=list(edge_counts[,1],edge_counts[,2]),sum,na.rm=TRUE)[,3],H,H)
diag(Edge) <- diag(Edge)/2

No_Edge <- matrix(aggregate(non_edge_counts[,3],by=list(non_edge_counts[,1],non_edge_counts[,2]),sum,na.rm=TRUE)[,3],H,H)
diag(No_Edge) <- diag(No_Edge)/2

a_n <- lowerTriangle(Edge,diag=TRUE)+a
b_bar_n <- lowerTriangle(No_Edge,diag=TRUE)+b

return(sum(lbeta(a_n,b_bar_n))-(H*(H+1)/2)*lbeta(a,b))
options(warn=0)

}

##########################################################################
## SIMULATION SCENARIOS FOR PEX-SBM WITH HPD, HNSP AND HDP HYPERPRIOR ####
##########################################################################

load("simulations.Rdata")

# setting the dimensions and other useful quantities
n_iter  <- 10000
burn_in <- 2000
V <- 80
N  <- c(30,30,15,5)

# HPD hyperparameters setting
theta  <- 0.5
theta0 <- 4

# HNSP hyperparameters setting
sigma  <- 0.2
sigma0 <- 0.8

# HPD hyperprior parameters setting
alpha   <- 5
lambda  <- 10
alpha0  <- 12
lambda0 <- 3

# check traceplots of beta-binomial likelihood for the first of the 10 replicated datasets
r <- 1	
Y1 <- Y1_replicates[,,r]
Y2 <- Y2_replicates[,,r]

#---------------------------------------------------------------
#---------------------------------------------------------------
# fit pEx-SBM(HDP) to scenario 1
set.seed(r)
model1       <- pex_sbm(Y1,n_iter,N,1,1,theta,theta0)
z_post1      <- model1[[1]]

ll_HDP_1 <- rep(0,n_iter)
for (s in 1:n_iter){
ll_HDP_1[s] <-  log_pY_z(Y1,z_post1[,s],1,1)	
if (s%%2500 == 0){print(paste("Iteration:", s))}
}

#---------------------------------------------------------------
# fit pEx-SBM(HDP) to scenario 2
set.seed(r)
model2       <- pex_sbm(Y2,n_iter,N,1,1,theta,theta0)
z_post2      <- model2[[1]]

ll_HDP_2 <- rep(0,n_iter)
for (s in 1:n_iter){
ll_HDP_2[s] <-  log_pY_z(Y2,z_post2[,s],1,1)	
if (s%%2500 == 0){print(paste("Iteration:", s))}
}

#---------------------------------------------------------------
#---------------------------------------------------------------
# fit pEx-SBM(HNSP) to scenario 1
set.seed(r)
model1       <- pex_sbm(Y1,n_iter,N,1,1,sigma=sigma, sigma0=sigma0, prior="HNSP")
z_post1      <- model1[[1]]

ll_HNSP_1 <- rep(0,n_iter)
for (s in 1:n_iter){
ll_HNSP_1[s] <-  log_pY_z(Y1,z_post1[,s],1,1)	
if (s%%2500 == 0){print(paste("Iteration:", s))}
}

#---------------------------------------------------------------
# fit pEx-SBM(HNSP) to scenario 2
set.seed(r)
model2       <- pex_sbm(Y2,n_iter,N,1,1,sigma=sigma, sigma0=sigma0, prior="HNSP")
z_post2      <- model2[[1]]

ll_HNSP_2 <- rep(0,n_iter)
for (s in 1:n_iter){
ll_HNSP_2[s] <-  log_pY_z(Y2,z_post2[,s],1,1)	
if (s%%2500 == 0){print(paste("Iteration:", s))}
}

#---------------------------------------------------------------
#---------------------------------------------------------------
# fit pEx-SBM(HDP hyperprior) to scenario 1
set.seed(r)
model1       <- hdpgamma_sbm(Y1,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post1      <- model1[[1]]

ll_HDP_HYP_1 <- rep(0,n_iter)
for (s in 1:n_iter){
ll_HDP_HYP_1[s] <-  log_pY_z(Y1,z_post1[,s],1,1)	
if (s%%2500 == 0){print(paste("Iteration:", s))}
}

#---------------------------------------------------------------
# fit pEx-SBM(HDP hyperprior) to scenario 2
set.seed(r)
model2       <- hdpgamma_sbm(Y2,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post2      <- model2[[1]]

ll_HDP_HYP_2 <- rep(0,n_iter)
for (s in 1:n_iter){
ll_HDP_HYP_2[s] <-  log_pY_z(Y2,z_post2[,s],1,1)	
if (s%%2500 == 0){print(paste("Iteration:", s))}
}

##########################################################################
################# SAVE BETA-BINOMIAL LOG-LIKELIHOODS #####################
##########################################################################

save(ll_HDP_1,ll_HDP_2,ll_HNSP_1,ll_HNSP_2,ll_HDP_HYP_1,ll_HDP_HYP_2,file="output_TRACEPLOTS.RData")

##########################################################################
################# FIGURES ################################################
##########################################################################
load("output_TRACEPLOTS.RData")

ll <- c(ll_HDP_1,ll_HDP_2,ll_HNSP_1,ll_HNSP_2,ll_HDP_HYP_1,ll_HDP_HYP_2)
group_1 <- c(rep("H-DP",n_iter*2),rep("H-NSP",n_iter*2),rep("H-DP hyp",n_iter*2))
group_2 <- c(rep("Scenario 1",n_iter),rep("Scenario 2",n_iter),rep("Scenario 1",n_iter),rep("Scenario 2",n_iter),rep("Scenario 1",n_iter),rep("Scenario 2",n_iter))

data    <- data.frame(ll)
data$g1 <- group_1
data$g1 <- factor(data$g1,levels=c("H-DP","H-NSP","H-DP hyp"))
data$g2 <- group_2
data$x <- rep(c(1:n_iter),6)

library(ggh4x)

df_scales <- data.frame(
  Panel = c("Scenario 1", "Scenario 2"),
  ymin = c(-1820, -2080),
  ymax = c(-1790, -2000),
  n = c(3, 3)
  )
df_scales <- split(df_scales, df_scales$Panel)

scales <- lapply(df_scales, function(x) {
  scale_y_continuous(limits = c(x$ymin, x$ymax), n.breaks = x$n)
})

ggplot(data=data, aes(x=x, y=ll)) + geom_line()+facet_grid(g2~ g1,scales="free")+
  ggh4x::facetted_pos_scales(   y = scales  )+  theme_bw()+ xlab("MCMC iterations")+ylab("Beta-Binomial log-likelihood")