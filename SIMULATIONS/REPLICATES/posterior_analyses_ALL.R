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
sourceCpp("JCDC.cpp")
source("casc.R")

##########################################################################
################# HPD ####################################################
##########################################################################
load("output_HDP.RData")

apply(waic,1,mean)

apply(vi,1,mean)

apply(expected_vi,1,mean)

apply(vi_bound,1,mean)

apply(H_estimated,1,mean)

apply(H_posterior,c(1,3),mean)

mean(misclass)

mean(pred_prob_err)

##########################################################################
################# HNSP ###################################################
##########################################################################
load("output_HNSP.RData")

apply(waic,1,mean)

apply(vi,1,mean)

apply(expected_vi,1,mean)

apply(vi_bound,1,mean)

apply(H_estimated,1,mean)

apply(H_posterior,c(1,3),mean)

mean(misclass)

mean(pred_prob_err)

##########################################################################
################# HPD Hyperprior #########################################
##########################################################################
load("output_HDP_hyperprior.RData")

apply(waic,1,mean)

apply(vi,1,mean)

apply(expected_vi,1,mean)

apply(vi_bound,1,mean)

apply(H_estimated,1,mean)

apply(H_posterior,c(1,3),mean)

mean(misclass)

mean(pred_prob_err)

##########################################################################
################# ESBM ###################################################
##########################################################################
load("output_ESBM.RData")

apply(waic,1,mean)

apply(vi,1,mean)

apply(expected_vi,1,mean)

apply(vi_bound,1,mean)

apply(H_estimated,1,mean)

apply(H_posterior,c(1,3),mean)

##########################################################################
################# COMPETITORS ############################################
##########################################################################
load("output_COMPETITORS.RData")

apply(vi,c(1,2),mean)

apply(H_estimated,c(1,2),mean)

##########################################################################
################# FIGURES ################################################
##########################################################################
load("simulations.Rdata")

# setting the dimensions and other useful quantities
V <- 80
N  <- c(30,30,15,5)

#---------------------------------------------------------------------------
# heatmap Y1 (first example)
Y <- Y1_replicates[,,7]

row_plot_Y <- as.data.frame(as.factor(z))
col_plot_Y <- as.data.frame(as.factor(lay))

names(row_plot_Y) <- "z"
names(col_plot_Y) <- "lay"
rownames(Y) <- rownames(row_plot_Y)
colnames(Y) <- rownames(col_plot_Y)

mycolors1 <- c(brewer.pal(9,"Reds")[6],brewer.pal(9,"Blues")[6],brewer.pal(9,"Greens")[6],brewer.pal(9,"Purples")[6])
names(mycolors1) <- unique(col_plot_Y$lay)

mycolors2 <- c(brewer.pal(9,"Greys")[1],brewer.pal(9,"Reds")[c(3,5)],brewer.pal(9,"Blues")[c(3,5)],brewer.pal(9,"Greens")[c(3,5)],brewer.pal(9,"Purples")[4])
names(mycolors2) <- unique(row_plot_Y$z)

mycolors <- list(z=mycolors2,lay=mycolors1)

Network1 <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_col = col_plot_Y,annotation_names_row=F,annotation_names_col=F, show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors = mycolors,gaps_row=c(which(diff(lay)!=0)),gaps_col=c(which(diff(lay)!=0)))

#---------------------------------------------------------------------------
# heatmap true edge probability matrix (scenario 1)
z_d <- as.matrix(dummy_cols(z)[,-1])
Probability_matrix <- z_d%*%Phi%*%t(z_d)
diag(Probability_matrix) <- 0

names(row_plot_Y) <- "z"
names(col_plot_Y) <- "lay"
rownames(Probability_matrix) <- rownames(row_plot_Y)
colnames(Probability_matrix) <- rownames(col_plot_Y)

mycolors1 <- c(brewer.pal(9,"Reds")[6],brewer.pal(9,"Blues")[6],brewer.pal(9,"Greens")[6],brewer.pal(9,"Purples")[6])
names(mycolors1) <- unique(col_plot_Y$lay)

mycolors2 <- c(brewer.pal(9,"Greys")[1],brewer.pal(9,"Reds")[c(3,5)],brewer.pal(9,"Blues")[c(3,5)],brewer.pal(9,"Greens")[c(3,5)],brewer.pal(9,"Purples")[4])
names(mycolors2) <- unique(row_plot_Y$z)

mycolors <- list(z=mycolors2,lay=mycolors1)

Network2 <- pheatmap(Probability_matrix,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),breaks=seq(0,1,by=1/30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_col = col_plot_Y,annotation_names_row=F,annotation_names_col=F, show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,
annotation_legend=F,annotation_colors = mycolors,gaps_row=c(which(diff(lay)!=0)),gaps_col=c(which(diff(lay)!=0)))

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
g <- grid.arrange(Network2[[4]],Network1[[4]],nrow=1,ncol=2,vp=viewport(width=1, height=1))
g_probal <- cowplot::ggdraw(g)+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[3]))
print(g_probal)
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# heatmap Y2 (first example)
Y <- Y2_replicates[,,7]

row_plot_Y <- as.data.frame(as.factor(z))
col_plot_Y <- as.data.frame(as.factor(lay))

names(row_plot_Y) <- "z"
names(col_plot_Y) <- "lay"
rownames(Y) <- rownames(row_plot_Y)
colnames(Y) <- rownames(col_plot_Y)

mycolors1 <- c(brewer.pal(9,"Reds")[6],brewer.pal(9,"Blues")[6],brewer.pal(9,"Greens")[6],brewer.pal(9,"Purples")[6])
names(mycolors1) <- unique(col_plot_Y$lay)

mycolors2 <- c(brewer.pal(9,"Greys")[1],brewer.pal(9,"Reds")[c(3,5)],brewer.pal(9,"Blues")[c(3,5)],brewer.pal(9,"Greens")[c(3,5)],brewer.pal(9,"Purples")[4])
names(mycolors2) <- unique(row_plot_Y$z)

mycolors <- list(z=mycolors2,lay=mycolors1)

Network3 <- pheatmap(Y,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_col = col_plot_Y,annotation_names_row=F,annotation_names_col=F, show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,annotation_legend=F,annotation_colors = mycolors,gaps_row=c(which(diff(lay)!=0)),gaps_col=c(which(diff(lay)!=0)))

#---------------------------------------------------------------------------
# heatmap true edge probability matrix (scenario 2)
z_d <- as.matrix(dummy_cols(z)[,-1])
Probability_matrix <- z_d%*%Phi_noise%*%t(z_d)
diag(Probability_matrix) <- 0

names(row_plot_Y) <- "z"
names(col_plot_Y) <- "lay"
rownames(Probability_matrix) <- rownames(row_plot_Y)
colnames(Probability_matrix) <- rownames(col_plot_Y)

mycolors1 <- c(brewer.pal(9,"Reds")[6],brewer.pal(9,"Blues")[6],brewer.pal(9,"Greens")[6],brewer.pal(9,"Purples")[6])
names(mycolors1) <- unique(col_plot_Y$lay)

mycolors2 <- c(brewer.pal(9,"Greys")[1],brewer.pal(9,"Reds")[c(3,5)],brewer.pal(9,"Blues")[c(3,5)],brewer.pal(9,"Greens")[c(3,5)],brewer.pal(9,"Purples")[4])
names(mycolors2) <- unique(row_plot_Y$z)

mycolors <- list(z=mycolors2,lay=mycolors1)

Network4 <- pheatmap(Probability_matrix,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(30),breaks=seq(0,1,by=1/30),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_col = col_plot_Y,annotation_names_row=F,annotation_names_col=F, show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,
annotation_legend=F,annotation_colors = mycolors,gaps_row=c(which(diff(lay)!=0)),gaps_col=c(which(diff(lay)!=0)))


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
g <- grid.arrange(Network4[[4]],Network3[[4]],nrow=1,ncol=2,vp=viewport(width=1, height=1))
g_probal <- cowplot::ggdraw(g)+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[3]))
print(g_probal)
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
