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
library(plyr)

# remove all except for functions
rm(list = setdiff(ls(), lsf.str()))

# load source functions
source("pex-sbm.R")
source("numcl.R")

##########################################################################
################# APPLICATION ############################################
##########################################################################

load("crime_net.RData")

# order wrt locali (i.e. layers) and set the dimensions and other useful quantities
V  <- dim(Y)[1]
df <- data.frame(Locale,index=c(1:V))
df <- arrange(df,Locale)

Y <- Y[df$index,df$index]
Role <- Role[df$index]
RoleLocale <- RoleLocale[df$index]

d <- nlevels(as.factor(df$Locale))
N <- as.numeric(table(df$Locale))
lay <- as.factor(df$Locale)
cumul <- cumsum(N)

n_iter  <- 10000
burn_in <- 2000

alpha   <- 10
lambda  <- 2.5
alpha0  <- 5
lambda0 <- 0.45
expnumcl_hdpgamma(N,alpha,lambda,alpha0,lambda0,1000)

theta_esbm <-  5
expnumcl_dp(V,theta_esbm)

###################################################################
# INFERENCE
###################################################################

#---------------------------------------------------------------
# fit pEx-SBM(HDP hyperprior)
set.seed(1)
model_HYP    <- hdpgamma_sbm(Y,n_iter,N,1,1,alpha,lambda,alpha0,lambda0)
z_post_HYP   <- model_HYP[[1]]
P_post_HYP   <- pr_cc(z_post_HYP)
z_hat_HYP    <- minVI(P_post_HYP,method="avg",max.k=20)$cl
credibleball(z_hat_HYP,t(z_post_HYP))[[5]]
table(z_hat_HYP)

LL <- matrix(nrow=V*(V-1)/2,ncol=n_iter-burn_in)
z_waic_HYP <- z_post_HYP[,(burn_in+1):n_iter]
for (t in 1:dim(z_waic_HYP)[2]){
  LL[,t]<-sampleLL(z_waic_HYP[,t],Y,1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
waic_HYP <- WAIC(LL)$WAIC
waic_HYP

#---------------------------------------------------------------
# fit ESBM
set.seed(1)
model_ESBM    <- esbm(Y,1,n_iter,"DP",alpha_PY = theta_esbm, sigma_PY = 0,x=as.numeric(lay), alpha_xi = rep(1,5))
z_post_ESBM   <- model_ESBM
P_post_ESBM   <- pr_cc(z_post_ESBM)
z_hat_ESBM    <- minVI(P_post_ESBM,method="avg",max.k=20)$cl
table(z_hat_ESBM)

LL <- matrix(nrow=V*(V-1)/2,ncol=n_iter-burn_in)
z_waic_ESBM <- z_post_ESBM[,(burn_in+1):n_iter]
for (t in 1:dim(z_waic_ESBM)[2]){
  LL[,t]<-sampleLL(z_waic_ESBM[,t],Y,1,1)
  if (t%%1000 == 0){print(paste("Iteration:", t))}
}
waic_ESBM <- WAIC(LL)$WAIC
waic_ESBM

###################################################################
# EDGE PREDICTION
###################################################################

# remove m nodes at random
m <- 10
Vrem <- V-m
set.seed(10)
rem <- randsample(V,m)
Yrem <- Y[-rem,-rem]
layrem <- as.numeric(lay[-rem])
remlay <- as.numeric(lay[rem])
Nrem <- as.numeric(table(layrem))

#---------------------------------------------------------------
# fit pEx-SBM(HDP)
set.seed(1)

model_rem <- hdpgamma_sbm(Yrem,n_iter,Nrem,1,1,alpha,lambda,alpha0,lambda0)
z_post <- model_rem[[1]]
w_post <- model_rem[[2]]
hpars  <- as.numeric(model_rem[[3]][2:3])
zhat <- minVI(pr_cc(z_post),method="avg",max.k=20)$cl
n_mc <- n_iter - burn_in

pred <- pred_conn_hyp(Yrem,layrem,remlay,model_rem,n_mc)
Ynew <- pred[[1]]
perm <- c(c(1:84)[-rem],rem)


Ycompare <- Y[perm,perm]
vec_pred<-c()
vec_data<-c()

for (i in 74:84){
	for (j in 1:(i-1)){
		vec_pred<-c(vec_pred,Ynew[i,j])
		vec_data<-c(vec_data,Ycompare[i,j])
	}
}

library(pROC)
roc_obj <- roc(vec_data, vec_pred)
AUC <- auc(roc_obj)
AUC


###################################################################
# GROUP PREDICTION
###################################################################

# remove m nodes at random
m <- 10
Vrem <- V-m
set.seed(10)
rem <- randsample(V,m)
Yrem <- Y[-rem,-rem]
layrem <- as.numeric(lay[-rem])
remlay <- as.numeric(lay[rem])
Nrem <- as.numeric(table(layrem))


# fit pEx-SBM(HDP hyperprior)
set.seed(1)
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

# Figure S.2
z_hat_prediction  <- minVI(P,method="avg",max.k=20)$cl
Locale_plot_pred<-c(df$Locale[-rem],df$Locale[rem])

sel_plot<-c(order(z_hat_prediction[1:74]),order(z_hat_prediction[75:84])+74)
P_plot<-P[sel_plot,sel_plot]
Locale_plot_pred<-Locale_plot_pred[sel_plot]

row_plot_Y <- as.data.frame(as.factor(Locale_plot_pred))
names(row_plot_Y) <- "Locale_plot_pred"
rownames(P_plot) <- rownames(row_plot_Y)

mycolors <- c(brewer.pal(10,"RdBu")[c(3,8)],brewer.pal(10,"PRGn")[c(8,3)],brewer.pal(9,"YlOrBr")[5])
names(mycolors) <- unique(row_plot_Y$Locale_plot_pred)

Predictive <- pheatmap(P_plot,color=colorRampPalette(brewer.pal(9,"Greys")[c(1,8)])(60),breaks=seq(0,1,by=1/60),cluster_cols = F, cluster_rows= F,annotation_row = row_plot_Y,annotation_names_row=F,annotation_names_col=F, show_rownames=F,show_colnames=F,legend=F,border_color=FALSE,
annotation_legend=F,annotation_colors = list(Locale_plot_pred=mycolors),gaps_row=c(74),gaps_col=c(74))


g <- grid.arrange(Predictive[[4]],nrow=1,ncol=1,vp=viewport(width=1, height=1))
g_probal <- cowplot::ggdraw(g)+ theme(plot.background =element_rect(fill=brewer.pal(9,"Greys")[3]))
print(g_probal)


###################################################################
# FIGURE 4
###################################################################

sel<-order(df$index)
z_crime<-z_hat_HYP
z_crime<-z_crime[sel]

load("crime_net.RData")

#load("Application/crime_net.RData")
V <- dim(Y)[1]

# note that Y must have diagonal equal to 0
diag(Y)

mycolors <- c(brewer.pal(10,"RdBu")[c(3,8)],brewer.pal(10,"PRGn")[c(8,3)],brewer.pal(9,"YlOrBr")[5],
brewer.pal(10,"RdBu")[c(3,8)],brewer.pal(10,"PRGn")[c(8,3)],brewer.pal(9,"YlOrBr")[5])

# transform the adjacency matrix into an igraph object
net_Y <- graph.adjacency(Y, mode=c("undirected"), weighted=NULL, diag=FALSE)

# compute the node betweenness to be used for the size of the nodes
betw <- betweenness(net_Y)
# node sizes are proportional to their betweenness 
# Note: for graphical purposes, we consider a monotone transformation of such a measure
V(net_Y)$size <- sqrt(betw/1.5+mean(betw))*1.2

# node colors indicate the presumed locale membership and role
V(net_Y)$color <- adjustcolor(mycolors[c(as.factor(RoleLocale))], alpha.f = .7)

# node shapes represent leadership role
V(net_Y)$shape <- c("circle","square")[c(as.factor(Role))]

# additional graphical settings
V(net_Y)$frame.color <- "black"
V(net_Y)$label <- "" 
E(net_Y)$color <- brewer.pal(9,"Greys")[4]

group_id<-list(which(z_crime==1),which(z_crime==2),which(z_crime==3),which(z_crime==4),which(z_crime==5),which(z_crime==6),which(z_crime==7),which(z_crime==8),which(z_crime==9),which(z_crime==10),which(z_crime==11),which(z_crime==12),which(z_crime==13),which(z_crime==14))

# plot Figure 4
plot(net_Y, rescale=F, layout=l*0.5,edge.curved=.3,edge.width=0.5, mark.groups = group_id,mark.col="#f0f0f0", mark.border="black")



