###########################################################################
# USEFUL FUNCTIONS TO ELICIT THE PRIOR ####################################
###########################################################################
### compile all the auxiliary functions in pex-sbm.R

library(gmp)

#---------------------------------------------------------------------------
### urn schemes

urn_DP <- function(tab,theta0){
  return(c(tab,theta0))
}

urn_NSP <- function(tab,sigma0){
  K <- length(tab)
  return(c(tab-sigma0,K*sigma0))
}

joint_HDP <- function(l, q, N, theta, theta0){
  P <- q/(theta+N-1)
  P <- rbind(P,rep(0,ncol(q)))
  newtab <- c(colSums(l),theta0)*theta/(theta+N-1)/(theta0+sum(l))
  P <- cbind(P,newtab)
  vecP<-as.vector(P)
  return(vecP)
}

joint_HNSP <- function(l, q, N, j, sigma, sigma0){
  qmin <- q-sigma
  qmin <- qmin*(qmin>0)
  P <- qmin/(N-1)
  P <- rbind(P,rep(0,ncol(q)))
  newtab <- c(colSums(l)-sigma0,nrow(q)*sigma0)*(sum(l[j,])*sigma)/(N-1)/(sum(l))
  P <- cbind(P,newtab)
  vecP<-as.vector(P)
  return(vecP)
}

#---------------------------------------------------------------------------
### functions to simulate z from the prior

#------------------
# DP
CRP <- function(n,prior,theta,sigma){
  z <- 1
  for(i in 2:n){
    
    if (prior=="DP"){
      urn <- function(tab){return(urn_DP(tab,theta))}
    } else if (prior=="NSP"){
      urn <- function(tab){return(urn_NSP(tab,sigma))}
    }
    
    p <- urn(as.vector(table(z)))
    znew <- which(rmultinom(1,1,p)==1)
    z <- c(z,znew)
  }
  return(z)
}

#------------------
# HDP/HNSP
CRF <- function(lay,prior,theta,theta0,sigma,sigma0){
  
  if (prior=="HDP"){
    joint<-function(l, q, N, j){return(joint_HDP(l, q, N, theta, theta0))}
    urn <- function(tab){return(urn_DP(tab,theta0))}
  } else if (prior=="HNSP"){
    joint<-function(l, q, N, j){return(joint_HNSP(l, q, N, j, sigma, sigma0))}
    urn <- function(tab){return(urn_NSP(tab,sigma0))}
  }
  
  n <- length(lay)
  
  z <- 1
  w <- 1
  
  #sample sequentially allocations for new nodes via the urn scheme
  for (i in 2:n){
    
    #get summaries for (z_1,...,z_i)
    out <- getsummary(z,w,lay[1:(i-1)])
    l <- out[[1]]
    q <- out[[2]]
    H <- max(z) 
    T <- length(unique(w[lay[1:(i-1)]==lay[i]])) #number of subgroups in layer lay[i]
    
    #if T=0, the i-th node is entering a new layer, so it creates a new subgroup and samples a cluster
    #from the simple urn scheme
    if(T==0){
      w <- c(w,1)
      p <- urn(as.vector(table(z)))
      znew <- which(rmultinom(1,1,p)==1)
      z <- c(z,znew)
    }else{
      #sample jointly allocation and suballocation for i-th node
      p <- joint(l,q[[lay[i]]],i,lay[i])
      new_sample <- matrix(rmultinom(1,1,p),H+1,T+1)
      wnew <- which(colSums(new_sample)==1)
      znew <- which(rowSums(new_sample)==1)
      
      #update z, w
      z <- c(z,znew)
      w <- c(w,wnew)}
  }
  
  return(z)
}

#------------------
# HDP with gamma hyperpriors
CRF_hp <- function(lay,alpha,lambda,alpha0,lambda0){
  
  joint<-function(l, q, N, j){return(joint_HDP(l, q, N, theta, theta0))}
  urn <- function(tab){return(urn_DP(tab,theta0))}
  
  theta  <- rgamma(1,alpha,lambda)
  theta0 <- rgamma(1,alpha0,lambda0)
  
  n <- length(lay)
  
  z <- 1
  w <- 1
  
  #sample sequentially allocations for new nodes via the urn scheme
  for (i in 2:n){
    
    #get summaries for (z_1,...,z_i)
    out <- getsummary(z,w,lay[1:(i-1)])
    l <- out[[1]]
    q <- out[[2]]
    H <- max(z) 
    T <- length(unique(w[lay[1:(i-1)]==lay[i]])) #number of subgroups in layer lay[i]
    
    #if T=0, the i-th node is entering a new layer, so it creates a new subgroup and samples a cluster
    #from the simple urn scheme
    if(T==0){
      w <- c(w,1)
      p <- urn(as.vector(table(z)))
      znew <- which(rmultinom(1,1,p)==1)
      z <- c(z,znew)
    }else{
      #sample jointly allocation and suballocation for i-th node
      p <- joint(l,q[[lay[i]]],i,lay[i])
      new_sample <- matrix(rmultinom(1,1,p),H+1,T+1)
      wnew <- which(colSums(new_sample)==1)
      znew <- which(rowSums(new_sample)==1)
      
      #update z, w 
      z <- c(z,znew)
      w <- c(w,wnew)}
  }
  
  return(z)
}

#---------------------------------------------------------------------------
### functions to simulate the number of non-empty clusters from the prior

#------------------
# DP

K_dp <- function(n,theta,mc_iter){
  K <- rep(0,mc_iter)
  for (i in 1:mc_iter){
    z <- CRP(n, prior="DP",theta)
    K[i] <- length(unique(z))
  }
  return(K)
}

#------------------
# HDP
K_hdp <- function(lay,theta,theta0,mc_iter){
  K <- rep(0,mc_iter)
  for (i in 1:mc_iter){
    laymix <- sample(lay)
    laymix <- as.numeric(fct_relabel(fct_inorder(as.factor(laymix)),~paste0(c(1:length(unique(laymix)))))) 
    z <- CRF(laymix, prior="HDP",theta,theta0)
    K[i] <- length(unique(z))
  }
  return(K)
}

#------------------
# HNSP
K_hnsp <- function(lay,sigma,sigma0,mc_iter){
  K <- rep(0,mc_iter)
  for (i in 1:mc_iter){
    laymix <- sample(lay)
    laymix <- as.numeric(fct_relabel(fct_inorder(as.factor(laymix)),~paste0(c(1:length(unique(laymix)))))) 
    z <- CRF(laymix, prior="HNSP",theta=NA,theta0=NA,sigma,sigma0)
    K[i] <- length(unique(z))
  }
  return(K)
}

#------------------
# HDP with gamma hyperpriors
K_hdpgamma <- function(lay,alpha,lambda,alpha0,lambda0,mc_iter){
  K <- rep(0,mc_iter)
  for (i in 1:mc_iter){
    laymix <- sample(lay)
    laymix <- as.numeric(fct_relabel(fct_inorder(as.factor(laymix)),~paste0(c(1:length(unique(laymix)))))) 
    z <- CRF_hp(laymix, alpha,lambda,alpha0,lambda0)
    K[i] <- length(unique(z))
  }
  return(K)
}


#---------------------------------------------------------------------------
### functions to study the prior on H and its expectation 
### with one marginalization step that reduces Monte Carlo variance

#------------------
# prob distr of the number of clusters with DP
pnumcl_dp <- function(nj,theta){

  out <-rep(0,nj)
  for(i in 1:nj){
    out[i] <- as.numeric(abs(Stirling1(nj,i)))*theta^i*gamma(theta)/gamma(theta+nj)
  }
  return(out)
}

#------------------
# expected number of clusters dp (closed-form)
expnumcl_dp <- function(n,theta){
  sum(theta/(theta+seq(1:n)-1))
}

#------------------
# generalized factorial (needed for NSP)
genfact<-function(n,sigma){
  out<-matrix(0,n,n)
  diag(out)<-sigma^(c(1:n))
  for(j in 1:(n-1)){
    out[j+1,1]<-(j-sigma)*out[j,1]}
  for(j in 2:(n-1)){
    for(k in 2:j){
      out[j+1,k]<-(j-k*sigma)*out[j,k]+sigma*out[j,k-1]}}
  return(out)
}

#------------------
# prob distr of the number of clusters with NSP
pnumcl_nsp <- function(V, sigma){
  factorial(c(0:(V-1)))*genfact(V,sigma)[V,]/factorial(V-1)/sigma
}

#------------------
# expected number of clusters HDP
expnumcl_hdp <- function(V,theta, theta0, n_mc){
  #V is the vector of layer numerosities
  d  <- length(V)
  mc <- 0
  probs  <- sapply(V,pnumcl_dp,theta,simplify=F)
  for (r in 1:n_mc){
    sample <- sapply(probs, function(iter) rmultinom(1,1, iter),simplify=F)
    ell    <- sum(sapply(sample, function(iter) which(iter==1)))
    mc     <- mc+sum(theta0/(theta0+c(1:ell)-1))
  }
  return(mc/n_mc)
}

#------------------
# expected number of clusters HNSP
expnumcl_hnsp <- function(V,sigma, sigma0, n_mc){
  #V is the vector of layer numerosities
  d  <- length(V)
  mc <- 0
  probs  <- sapply(V,pnumcl_nsp,sigma,simplify = F)
  for (r in 1:n_mc){
    
    sample <- sapply(probs, function(iter) rmultinom(1,1, iter),simplify=F)
    ell    <- sum(sapply(sample, function(iter) which(iter==1)))
    mc     <- mc+sum(c(1:ell)*pnumcl_nsp(ell,sigma0))
  }
  return(mc/n_mc)
}

#------------------
# expected number of clusters HDP with hyperprior
expnumcl_hdpgamma <- function(V,alpha, lambda, alpha0, lambda0, n_mc){
  #V is the vector of layer numerosities
  d  <- length(V)
  mc <- 0
  for (r in 1:n_mc){
    theta  <- rgamma(1,alpha,rate=lambda)
    probs  <- sapply(V,pnumcl_dp,theta,simplify=F)
    sample <- sapply(probs, function(iter) rmultinom(1,1, iter),simplify=F)
    ell    <- sum(sapply(sample, function(iter) which(iter==1)))
    theta0 <- rgamma(1,alpha0,rate=lambda0)
    mc     <- mc+sum(theta0/(theta0+c(1:ell)-1))
  }
  return(mc/n_mc)
}

#------------------
# distribution over number of clusters HNSP
numcl_hnsp <- function(V,sigma,sigma0,n_mc){
  #V is the vector of layer numerosities
  d  <- length(V)
  mc <- rep(0,sum(V))
  for (r in 1:n_mc){
    probs  <- sapply(V,pnumcl_nsp,sigma,simplify=F)
    sample <- sapply(probs, function(iter) rmultinom(1,1, iter),simplify=F)
    ell    <- sum(sapply(sample, function(iter) which(iter==1)))
    mc     <- mc+c(pnumcl_nsp(ell,sigma0),rep(0,sum(V)-ell))
  }
  return(mc/n_mc)
  
}

#------------------
# distribution over number of clusters HDP
numcl_hdp <- function(V,theta,theta0,n_mc){
  #V is the vector of layer numerosities
  d  <- length(V)
  mc <- rep(0,sum(V))
  for (r in 1:n_mc){
    probs  <- sapply(V,pnumcl_dp,theta)
    sample <- sapply(probs, function(iter) rmultinom(1,1, iter))
    ell    <- sum(sapply(sample, function(iter) which(iter==1)))
    mc     <- mc+c(pnumcl_dp(ell,theta0),rep(0,sum(V)-ell))
  }
  return(mc/n_mc)
  
}

#------------------
# distribution over number of clusters HDP hyperprior
numcl_hdpgamma <- function(V,alpha, lambda, alpha0, lambda0,n_mc){
  #V is the vector of layer numerosities
  d  <- length(V)
  mc <- rep(0,sum(V))
  for (r in 1:n_mc){
  	theta  <- rgamma(1,alpha,rate=lambda)
    probs  <- sapply(V,pnumcl_dp,theta)
    sample <- sapply(probs, function(iter) rmultinom(1,1, iter))
    ell    <- sum(sapply(sample, function(iter) which(iter==1)))
    theta0 <- rgamma(1,alpha0,rate=lambda0)
    mc     <- mc+c(pnumcl_dp(ell,theta0),rep(0,sum(V)-ell))
  }
  return(mc/n_mc)
  
}