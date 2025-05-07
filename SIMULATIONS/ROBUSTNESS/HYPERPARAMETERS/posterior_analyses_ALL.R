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

##########################################################################
################# HDP ####################################################
##########################################################################
load("output_HDP.RData")

apply(expected_vi,1,mean)

##########################################################################
################# HNSP ###################################################
##########################################################################
load("output_HNSP.RData")

apply(expected_vi,1,mean)

##########################################################################
################# HPD Hyperprior #########################################
##########################################################################
load("output_HDP_hyperprior.RData")

apply(expected_vi,1,mean)

##########################################################################
################# ESBM ###################################################
##########################################################################
load("output_ESBM.RData")

apply(expected_vi,1,mean)