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

load("output_HDP.RData")
dati_HDP<-melt(vi)[,-2]
dati_HDP$method<-rep("H-DP",140)

load("output_HNSP.RData")
dati_HNSP<-melt(vi)[,-2]
dati_HNSP$method<-rep("H-NSP",140)

load("output_HDP_hyperprior.RData")
dati_HDP_hyperprior<-melt(vi)[,-2]
dati_HDP_hyperprior$method<-rep("H-DP hyp",140)

load("output_CASC.RData")
dati_CASC<-melt(vi)[,-2]
dati_CASC$method<-rep("CASC",140)

dati_plot<-rbind(dati_HDP,dati_HDP_hyperprior,dati_HNSP,dati_CASC)

dati_plot$X1<-factor(dati_plot$X1)
dati_plot$X3<-factor(dati_plot$X3,levels=c("7","6","5","4","3","2","1"))

ggplot(dati_plot, aes(x = X3, y = value,fill=method))+geom_boxplot(outlier.size=0.5) + facet_wrap(~ X1,scales = "free_x",labeller=as_labeller(c(`1` = "Scenario 1", `2`="Scenario 2")))+
  theme_bw()+ scale_x_discrete(labels=c("7" = "40", "6" = "60","5" = "80", "4"="100","3"="120","2"="140","1"="160" )) + ylab(expression(VI(hat("z"),z[0])))+scale_fill_manual(values = c("#ffffff","#e6e6e6","#bfbfbf","#999999"))+ xlab("Number of nodes")+  theme(legend.title=element_blank())