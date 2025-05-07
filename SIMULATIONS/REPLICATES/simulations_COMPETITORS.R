#  load useful libraries
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
sourceCpp("JCDC.cpp")
source("casc.R")

##########################################################################
################# SIMULATION SCENARIOS FOR COMPETITORS ###################
##########################################################################

load("simulations.Rdata")

# setting the dimensions and other useful quantities
V <- 80
N  <- c(30,30,15,5)

# for JCDC
similarity <- matrix(NA,V,V)
for (i in 1:V)
  for (j in 1:V)
    similarity[i,j]<-as.integer(lay[i]==lay[j])

# for CASC
X_casc <- as.matrix(dummy_cols(lay)[,-1])

# define the matrices containing the performance measures
vi            <- array(0,c(2,5,10))
z_estimated   <- array(0,c(2,5,10,80))
H_estimated   <- array(0,c(2,5,10))


for (r in 1:10){
		
Y1 <- Y1_replicates[,,r]
Y2 <- Y2_replicates[,,r]



#---------------------------------------------------------------
# Louvain scenario 1
set.seed(r)
net <- graph.adjacency(Y1, mode=c("undirected"), weighted=NULL, diag=FALSE)
Louv <- cluster_louvain(net)$membership

vi[1,1,r] <- VI(z,t(Louv))
z_estimated[1,1,r,] <- Louv
H_estimated[1,1,r] <- max(Louv)

#---------------------------------------------------------------
# Louvain scenario 2
set.seed(r)
net <- graph.adjacency(Y2, mode=c("undirected"), weighted=NULL, diag=FALSE)
Louv <- cluster_louvain(net)$membership

vi[2,1,r] <- VI(z,t(Louv))
z_estimated[2,1,r,] <- Louv
H_estimated[2,1,r] <- max(Louv)



#---------------------------------------------------------------
# SBM scenario 1
set.seed(r)
greed_out <- greed(Y1,K=8,model=new("Sbm",alpha=1,type="undirected"),verbose=FALSE)
g_sbm <- greed_out@cl

vi[1,2,r] <- VI(z,t(g_sbm))
z_estimated[1,2,r,] <- g_sbm
H_estimated[1,2,r] <- max(g_sbm)

#---------------------------------------------------------------
# SBM scenario 2
set.seed(r)
greed_out <- greed(Y2,K=8,model=new("Sbm",alpha=1,type="undirected"),verbose=FALSE)
g_sbm <- greed_out@cl

vi[2,2,r] <- VI(z,t(g_sbm))
z_estimated[2,2,r,] <- g_sbm
H_estimated[2,2,r] <- max(g_sbm)



#---------------------------------------------------------------
# JCDC scenario 1
A <- Y1
K <- 8

# initialization
D.inv = diag(1./(sqrt(apply(A, 1, sum))+1e-7));
Laplacian = D.inv %*% A %*% D.inv;
L.svd = svd(Laplacian);
U.K = L.svd$u[, 1:K];
spec.cluster = kmeans(U.K, K, nstart=10)$cluster;

G.fit = array(0, c(V, K));
for(k in 1:K){
  G.fit[spec.cluster==k, k] = 1;
}

#--------------------
# tuning parameter 5
set.seed(r)
W_max <- 5
p = dim(similarity)[3];	if(length(dim(similarity))==2){p = 1;}

result <- JCDC(A, similarity, p, G.fit, 1, K, 20, 30, 20, W_max, 2, 0, 1);
cluster_JCDC = as.vector(result$G.fit %*% (1:K));
jcdc_5 <- cluster_JCDC

vi[1,3,r] <- VI(z,t(jcdc_5))
z_estimated[1,3,r,] <- jcdc_5
H_estimated[1,3,r] <- max(jcdc_5)

#--------------------
# tuning parameter 1.5
set.seed(r)
W_max <- 1.5
p = dim(similarity)[3];	if(length(dim(similarity))==2){p = 1;}

result <- JCDC(A, similarity, p, G.fit, 1, K, 20, 30, 20, W_max, 2, 0, 1);
cluster_JCDC = as.vector(result$G.fit %*% (1:K));
jcdc_1_5 <- cluster_JCDC

vi[1,4,r] <- VI(z,t(jcdc_1_5))
z_estimated[1,4,r,] <- jcdc_1_5
H_estimated[1,4,r] <- max(jcdc_1_5)

#---------------------------------------------------------------
# JCDC scenario 2
A <- Y2
K <- 8

# initialization
D.inv = diag(1./(sqrt(apply(A, 1, sum))+1e-7));
Laplacian = D.inv %*% A %*% D.inv;
L.svd = svd(Laplacian);
U.K = L.svd$u[, 1:K];
spec.cluster = kmeans(U.K, K, nstart=10)$cluster;

G.fit = array(0, c(V, K));
for(k in 1:K){
  G.fit[spec.cluster==k, k] = 1;
}

#--------------------
# tuning parameter 5
set.seed(r)
W_max <- 5
p = dim(similarity)[3];	if(length(dim(similarity))==2){p = 1;}

result <- JCDC(A, similarity, p, G.fit, 1, K, 20, 30, 20, W_max, 2, 0, 1);
cluster_JCDC = as.vector(result$G.fit %*% (1:K));
jcdc_5 <- cluster_JCDC

vi[2,3,r] <- VI(z,t(jcdc_5))
z_estimated[2,3,r,] <- jcdc_5
H_estimated[2,3,r] <- max(jcdc_5)

#--------------------
# tuning parameter 1.5
set.seed(r)
W_max <- 1.5
p = dim(similarity)[3];	if(length(dim(similarity))==2){p = 1;}

result <- JCDC(A, similarity, p, G.fit, 1, K, 20, 30, 20, W_max, 2, 0, 1);
cluster_JCDC = as.vector(result$G.fit %*% (1:K));
jcdc_1_5 <- cluster_JCDC

vi[2,4,r] <- VI(z,t(jcdc_1_5))
z_estimated[2,4,r,] <- jcdc_1_5
H_estimated[2,4,r] <- max(jcdc_1_5)



#---------------------------------------------------------------
# CASC scenario 1 - Binkiewicz et al. (2017)
K <- 8
set.seed(r)
cluster_casc <- casc(Y1,X_casc,K)$cluster

vi[1,5,r] <- VI(z,t(cluster_casc))
z_estimated[1,5,r,] <- cluster_casc
H_estimated[1,5,r] <- max(cluster_casc)

#---------------------------------------------------------------
# CASC scenario 2 - Binkiewicz et al. (2017)
K <- 8
set.seed(r)
cluster_casc <- casc(Y2,X_casc,K)$cluster

vi[2,5,r] <- VI(z,t(cluster_casc))
z_estimated[2,5,r,] <- cluster_casc
H_estimated[2,5,r] <- max(cluster_casc)


}

save(vi,z_estimated,H_estimated, file="output_COMPETITORS.RData")
