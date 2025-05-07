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
################# SIMULATION SCENARIOS FOR CASC  #########################
##########################################################################

load("simulations_consistency.Rdata")

# define the matrices containing the performance measures
vi            <- array(0,c(2,10,7))

for (r in 1:7){
		
Y1_temp <- Y1_consistency[[r]]
Y2_temp <- Y2_consistency[[r]]

lay <- lay_consistency[[r]]
z <- z_consistency[[r]]
N <- c(table(lay))
X_casc <- as.matrix(dummy_cols(lay)[,-1])

for (s in 1:10){
	
Y1 <- Y1_temp[,,s]
Y2 <- Y2_temp[,,s]
V <- dim(Y1)[1]

#---------------------------------------------------------------
# CASC scenario 1 - Binkiewicz et al. (2017)
K <- 8
set.seed(s)
cluster_casc <- casc(Y1,X_casc,K)$cluster

vi[1,s,r] <- VI(z,t(cluster_casc))

#---------------------------------------------------------------
# CASC scenario 2 - Binkiewicz et al. (2017)
K <- 8
set.seed(s)
cluster_casc <- casc(Y2,X_casc,K)$cluster

vi[2,s,r] <- VI(z,t(cluster_casc))

print(s)}
print(r)}

save(vi,file="output_CASC.RData")
