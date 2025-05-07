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
################# ROBUSTNESS #############################################
##########################################################################
load("output_ROBUSTNESS.RData")

apply(vi,1,mean)
