library(reshape)
library(gdata)
library(igraph)
library(mcclust.ext)
library(RColorBrewer)
library(pheatmap)
library(gridExtra)
library(grid)
library(cowplot)
library(ggplot2)
library(coda)
library(dummy)
library(randnet)
library(greed)
library(LaplacesDemon)
library(tidyverse)
library(gmp)

####################################################################
################# MAIN FUNCTIONS FOR PEX-SBM #######################
####################################################################

#MCMC FOR SAMPLING FROM POSTERIOR DISTRIBUTION OF THE PARTITION

pex_sbm <- function(Y, n_iter, N, a, b, theta=NA, theta0=NA, sigma=NA, sigma0=NA, prior="HDP", z_init=c(1:nrow(Y))){
  #input: Y adjacency matrix
  #       n_iter number of mcmc iterations
  #       N vector of numerosities of the layers
  #       a,b beta pars
  #       theta, theta0, sigma, sigma0 prior hyperpars
  #       prior "HDP" or "HNSP"
  #       z_init initial allocation
  #output: a list of posterior allocation sample
  #                  posterior sub-allocation sample
  #                  vector of prior type and hyperpars
  #                  
  #WARNING: needs more than 1 node in each layer and more than 1 layer
  
  #check numerosities
  if (sum(N)!=nrow(Y)){
    stop("Number of nodes across layers does not agree with number of data rows")
  }
  
  #assign joint urn scheme according to prior selection
  if (prior=="HDP"){
    joint<-function(l, q, N, j){return(joint_HDP(l, q, N, theta, theta0))}
  } else if (prior=="HNSP"){
    joint<-function(l, q, N, j){return(joint_HNSP(l, q, N, j, sigma, sigma0))}
  }
  
  #total number of nodes
  V <- sum(N)
  #total number of layers
  d <- length(N)
  #layer indicator
  lay <- numeric(0)
  for(j in 1:d){
    lay <- c(lay,rep(j,N[j]))
  }
  
  # ----------------------------------------------
  # Initialization 
  # ----------------------------------------------
  
  #initialize matrix version of the allocation vector: Z[v,h]=1 if node v is in cluster h, Z[v,h]=0 otherwise
  Z <- vec2mat(z_init)
  
  #allocation and suballocation, each column stores an mcmc output 
  z_post <- matrix(NA,V,n_iter)
  t_post <- matrix(NA,V,n_iter)
  
  #initialize list version of suballocation vector: all singletons, coherent with every initialization of the allocations
  T <- vector("list", d)
  for (j in 1:d){
    T[[j]]<-c(1:N[j])
  }
  
  #initialize summaries l and q:
  # matrix l[j,h] = number of subgroups in layer j and cluster h
  # list q[[j]][h,t] = number of nodes in layer j in cluster h and subgroup t (t runs layer-wise)
  l<-matrix(0,d,ncol(Z))
  q <- vector("list", d)
  for (j in 1:d){
    l[j,]  <- colSums(as.matrix(Z[lay==j,]))
    q[[j]] <- t(as.matrix(Z[lay==j,]))
  }  
  
  
  #initialize count matrix of block connections
  temp   <- Y%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  
  # ----------------------------------------------
  # Gibbs sampler
  # ----------------------------------------------
  
  for (s in 1:n_iter){
    i<-0
    j<-1
    
    #running across nodes
    for (v in 1:V){
      
      #set layer and in-layer node counters
      i<-i+1
      if(i>N[j]){
        i<-1
        j<-j+1}
      
      ## Remove v from Z and l                                                  ##
      ## if the cluster containing node v has no other node, discard it as well ##
      
      #non-empty clusters after removing v
      nonempty_v <- which(colSums(as.matrix(Z[-v,])) > 0)
      
      if(ncol(Z) > 1){
        
        #reduce Z and l
        Z <- Z[, nonempty_v]
        l <- l[, nonempty_v]
        
        #if it's just one cluster, keep Z and l as matrices
        if (length(nonempty_v)==1){
          Z <- matrix(Z,V,1)
          l <- matrix(l,d,1)}
        
        #reduce the dimensions of the m_full matrix if we discard a cluster
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      
      # Z_v = allocation matrix excluding node v
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y[-v,v])
      
      # h_v = cluster of the removed node, empty if it was a singleton
      h_v         <- which(Z[v,] > 0)
      
      ## Remove node v from q and T                                                ##
      ## if the subgroup containing node v has no other node, discard it as well ##
      
      #matrix version of the suballocation vector of layer j
      Tmat<-vec2mat(T[[j]])
      #non-empty subgroups after removing v
      nonempty_vt <- which(colSums(as.matrix(Tmat[-i,]))>0)
      
      #reduce T and q[[j]]
      if(ncol(Tmat)>1){
        Tmat        <- as.matrix(Tmat[,nonempty_vt])
        q[[j]]      <- as.matrix(q[[j]][nonempty_v,nonempty_vt])
        if(length(nonempty_v)==1){q[[j]] <- t(q[[j]])}
      }
      
      #if the cluster was removed, reduce the rest of q accordingly
      if(length(h_v)!=1){
        for(k in 1:d){
          if(k!=j){
            q[[k]] <- as.matrix(q[[k]][nonempty_v,])
            if(length(nonempty_v)==1)
            {q[[k]] <- t(q[[k]])}
          }
        }
        
      }
      
      ## Compute the count matrices of block connections (and non-connections) after removing v ##
      
      # subgroup of the removed node v
      t_v <- which(Tmat[i,] > 0)
      
      if(length(h_v)==1){# if didn't close a cluster
        resid1       <- matrix(0,H,H)
        resid1[,h_v] <- r_v; resid1[h_v,] <- r_v #residual connections are just r_v
        m            <- m_full - resid1 # compute the m matrix by difference
        if(length(t_v)==1){# if didn't close a subgroup
          q[[j]][h_v,T[[j]][i]] <- q[[j]][h_v,T[[j]][i]]-1 #one less node in the subgroup
        }else{l[j,h_v] <- l[j,h_v]-1} #if closed a subgroup, one less subgroup in the cluster
      }else {m <- m_full} # if closed a cluster (and hence a subgroup), already deleted rows and columns from q, l and m
      
      #rewrite the subgroup vector allocation T from Tmat
      T[[j]] <- Tmat %*% c(1:ncol(Tmat))
      
      #number of subgroups excluding v
      ntab <- ncol(q[[j]])
      
      
      # m_bar = number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      
      #matrix versions of v_minus and r_v
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      
      ## Compute the probabilities and sample cluster and subgroup for v##
      
      #loglikelihood
      log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      
      #summing predictive probabilities from joint urn scheme
      log_p <- log(joint(l,q[[j]],N[j],j)) + rep(c(log_lhds_old, log_lhd_new),ntab+1)
      p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      #sampling allocation and suballocation for v
      new_sample <- matrix(rmultinom(1,1,p),H+1,ntab+1)
      t_new <- which(colSums(new_sample)==1)
      h_new <- which(rowSums(new_sample)==1)
      
      ## Adjust Z, H, r_v, m, T, q and l ##
      
      #write the subllocation in T
      T[[j]][i]<-t_new
      
      #if the old cluster of v still exists, delete the allocation of v in it from Z
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      
      #if new subgroup
      if(t_new > max(T[[j]][-i])){
        q[[j]] <- cbind(q[[j]],rep(0,H))
        #... and new cluster
        if(h_new == H+1){
          Z               <- cbind(Z,rep(0,V)) #add a column to Z
          Z[v,h_new] <- 1 #allocate v in the new column
          m               <- rbind(cbind(m,0),0) #add a column to m
          r_v             <- crossprod(Z[-v,], Y[-v,v]) #redefine the residual connections
          l               <- cbind(l,rep(0,d)) #add a column to l
          l[j,h_new] <- 1 #count a new subgroup
          for(k in 1:d){
            if(k==j){
              q[[j]] <- rbind(q[[j]],c(rep(0,ntab),1)) #add a column to q[[j]] and count a new node
            }else{q[[k]] <- rbind(q[[k]],rep(0,ncol(q[[k]])))} #add a column to the rest of q
          }
          H <- H + 1 #add a cluster
        }else{#if old cluster
          l[j,h_new]            <- l[j,h_new]+1 #add a new subgroup
          q[[j]][h_new,ntab+1]  <- 1 #count a new node in the new subgroup and old cluster
          Z[v, h_new]           <- 1} #allocate v
      }else{#if old subgroup
        Z[v, h_new]    <- 1  #allocate v
        q[[j]][h_new,t_new] <- q[[j]][h_new,t_new]+1} #add a new node in the old subgroup and old cluster
      
      # update m_full
      resid2              <- matrix(0,H,H)
      resid2[,h_new] <- r_v; resid2[h_new,] <- r_v
      m_full              <- m + resid2
    }
    #end running across nodes
    
    
    #store reordered allocation
    z          <- Z %*% c(1:ncol(Z))
    z_post[,s] <- as.numeric(fct_relabel(fct_inorder(as.factor(z)),~paste0(c(1:length(unique(z))))))
    
    
    #store reordered suballocation
    t <- T
    for(r in 1:d){
      t[[r]] <- as.numeric(fct_relabel(fct_inorder(as.factor(t[[r]])),~paste0(c(1:length(unique(t[[r]]))))))
    }
    t_post[,s] <- list2vec(t,V)
    
    #iteration counter
    if (s%%1000 == 0){print(paste("Iteration:", s))}
  }
  
  # ----------------------------------------------
  # End of Gibbs sampler
  # ----------------------------------------------
  
  #store model specifications and posterior samples in the output list
  prior <- c(prior, theta, theta0, sigma, sigma0)
  pars<-c(a,b)
  model<-list(z_post,t_post,prior,pars)
  
  return(model)
}


#JOINT URN SCHEMES FOR HDP AND HNSP

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

#CO-CLUSTERING PROBABILITIES FOR ALLOCATION PREDICTION

coclust_new <- function(l, q, jr, js, theta, theta0, sigma=NA, sigma0=NA, prior="HDP", lay){
  #input:  l,q summaries
  #        jr, js layer position of involved new nodes
  #        theta, theta0, sigma, sigma0  hyperparameters
  #        prior "HDP" or "HNSP"
  #        lay vector of counts of layers of old nodes
  #output: co-clustering probability between two new nodes
  lsum  <- sum(l)
  lrows <- colSums(l)
  lcols <- rowSums(l)
  H     <- dim(l)[2]
  if(jr==js){
    j   <- jr
    nj  <- rowSums(q[[j]])
    if(prior=="HDP"){
    out <- (sum(nj*(nj+1))+theta*(1+2/(theta0+lsum)*sum(nj*lrows))+
              theta^2/(theta0+lsum)/(theta0+lsum+1)*(sum(lrows*(lrows+1))+theta0))/(theta+lay[j])/(theta+lay[j]+1)
    }
    else{
    out <- (sum((nj-l[j,]*sigma)*(nj-l[j,]*sigma+1))+lcols[j]*sigma*((1-sigma)+2/lsum*sum((nj-l[j,]*sigma)*(lrows-sigma0)))+
            lcols[j]*(lcols[j]+1)*sigma^2/lsum/(lsum+1)*(sum((lrows-sigma0)*(lrows-sigma0+1))+H*sigma0*(1-sigma0)))/lay[j]/(lay[j]+1)
    }
  }
  else{
    njr <- rowSums(q[[jr]])
    njs <- rowSums(q[[js]])
    if(prior=="HDP"){
    out <- (sum(njr*njs)+theta/(theta0+lsum)*sum((njr+njs)*lrows)+
              theta^2/(theta0+lsum)/(theta0+lsum+1)*(sum(lrows*(lrows+1))+theta0))/
      (theta+lay[jr])/(theta+lay[js])}
    else{
    out <- (sum((njr-l[jr,]*sigma)*(njs-l[js,]*sigma))+
              sigma/lsum*(lcols[jr]*sum((njs-l[js,]*sigma)*(lrows-sigma0))+lcols[js]*sum((njr-l[jr,]*sigma)*(lrows-sigma0)))+
           lcols[jr]*lcols[js]*sigma^2/lsum/(lsum+1)*(sum((lrows-sigma0)*(lrows-sigma0+1)+H*sigma0*(1-sigma0))))/(lay[jr]-1)/(lay[js]-1)
    }
  }
  return(out)  
}

coclust_old <- function(l, q, j, z, theta, theta0, sigma=NA, sigma0=NA, prior="HDP", lay){
  #input:  l,q summaries
  #        j layer position of the new node
  #        z cluster allocation of the old node
  #        theta, theta0, sigma, sigma0  hyperparameters
  #        prior "HDP" or "HNSP"
  #        lay vector of counts of layers of old nodes
  #output: co-clustering probability between new node and old node
  lsum  <- sum(l)
  lrows <- colSums(l)
  lcols <- rowSums(l)
  nj  <- rowSums(q[[j]])
  if(prior=="HDP"){
    out <- lrows[z]/(theta0+lsum)*theta/(theta+lay[j])+nj[z]/(theta+lay[j])}
  else{
    out <- (lrows[z]-sigma0)/lsum*lcols[j]*sigma/(lay[j])+(nj[z]-l[j,z]*sigma)/lay[j]
  }
  
  return(out)  
}

#EDGE PREDICTION

pred_conn <- function(Y,lay,newlay,model,nmc){
  #input: Y observed adjacency matrix
  #       lay vector of layer allocations of nodes in Y
  #       newlay vector of layer allocations of new nodes
  #       model output of pex-sbm function
  #       nmc number of MC iterations
  #output: list(adjacency matrix completed with posterior edge probabilities of new nodes, vector of number of new discoveries)
  
  #unpack model characteristics and posterior allocation samples
  N <- dim(Y)[1]
  m <- length(newlay)
  zpost <- model[[1]]
  wpost <- model[[2]]
  prior <- model[[3]][1]
  ab <- model[[4]]
  hyperpars <- as.numeric(model[[3]][2:5])
  nmcmc <- dim(zpost)[[2]]
  
  #check number of mc iterations
  if(nmc>nmcmc){print(paste("Number of MC iterations too large! Using instead number of MCMC iterations:", nmcmc))
    nmc<-nmcmc}
  
  if (prior=="HDP"){
    joint<-function(l, q, N, j){return(joint_HDP(l, q, N, hyperpars[1], hyperpars[2]))}
  } else if (prior=="HNSP"){
    joint<-function(l, q, N, j){return(joint_HNSP(l, q, N, j, hyperpars[3], hyperpars[4]))}
  }
  
  #initialize vector of total and old number of clusters 
  H <- rep(0,nmc)
  Hold <- rep(0,nmc)
  
  #initialize extended adjacency matrix
  Ynew <- cbind(nmc*Y,matrix(0,N,m))
  Ynew <- rbind(Ynew,matrix(0,m,N+m))
  
  #obtain connection probabilities of new nodes for each MC iteration and sum them up
  for (n in 1:nmc){
    
    #consider the (burnin+n)-th joint posterior allocation sample
    z <- zpost[,nmcmc-nmc+n]
    w <- wpost[,nmcmc-nmc+n]
    laytemp <- lay
    
    #sample sequentially allocations for new nodes via the urn scheme
    for (i in 1:m){
      
      #get summaries for (z_1,...,z_N,...,z_(N+i-1))
      out <- getsummary(z,w,laytemp)
      l <- out[[1]]
      q <- out[[2]]
      ell <- newlay[i] #(N+i)-th node's layer
      Htemp <- max(z) 
      T <- max(w[laytemp==ell]) #number of subgroups in layer ell
      
      #sample jointly allocation and suballocation for (N+i)-th node
      p <- joint(l,q[[ell]],N,ell)
      new_sample <- matrix(rmultinom(1,1,p),Htemp+1,T+1)
      wnew <- which(colSums(new_sample)==1)
      znew <- which(rowSums(new_sample)==1)
      
      #update z, w and lay: (N+i)-th node is old now!
      z <- c(z,znew)
      w <- c(w,wnew)
      laytemp <- c(laytemp,newlay[i])
    }
    
    
    #get the matrices M0 and M0bar of counts of edges and non-edges between clusters
    zold <- z[1:N]
    Z    <- vec2mat(zold)
    temp <- Y%*%Z
    M0   <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
    Mtot <- as.matrix(as.numeric(table(zold)))%*%t(as.matrix(as.numeric(table(zold))))
    diag(Mtot) <- as.numeric(table(zold))*(as.numeric(table(zold))-1)/2
    M0bar <- Mtot - M0
    
    H[n] <- max(z)
    Hold[n] <- max(z[1:N])
    M0   <- cbind(M0,matrix(0,Hold[n],H[n]-Hold[n]))
    M0   <- rbind(M0,matrix(0,H[n]-Hold[n],H[n]))
    M0bar   <- cbind(M0bar,matrix(0,Hold[n],H[n]-Hold[n]))
    M0bar   <- rbind(M0bar,matrix(0,H[n]-Hold[n],H[n]))
   
    
    #compute the posterior beta expectation for each extended MC sample 
    for (i in 1:m){
      for(j in 1:(N+i-1)){
        Ynew[N+i,j]<-Ynew[N+i,j]+(ab[1]+M0[z[N+i],z[j]])/(ab[1]+ab[2]+M0[z[N+i],z[j]]+M0bar[z[N+i],z[j]])
        Ynew[j,N+i]<-Ynew[N+i,j]
      }
    }
    if(n%%1000==0){print(paste("Iteration:",n))}
  } 
  return(list(Ynew/nmc,H-Hold))
  
}


###HDP WITH GAMMA HYPERPRIORS###

hdpgamma_sbm <- function(Y, n_iter, N, a, b, alpha, lambda, alpha0, lambda0, z_init=c(1:nrow(Y))){
  #input: Y adjacency matrix
  #       n_iter number of mcmc iterations
  #       N vector of numerosities of the layers
  #       a,b beta pars
  #       alpha, lambda (rate), alpha0, lambda0 (rate) gamma hyperprior pars
  #       z_init initial allocation
  #output: a list of posterior allocation sample
  #                  posterior sub-allocation sample
  #                  posterior of hyperparameters
  #                  vector of prior type and hyperpars
  #                  
  #WARNING: needs more than 1 node in each layer and more than 1 layer
  
  #check numerosities
  if (sum(N)!=nrow(Y)){
    stop("Number of nodes across layers does not agree with number of data rows")
  }
  
  #total number of nodes
  V <- sum(N)
  #total number of layers
  d <- length(N)
  #layer indicator
  lay <- numeric(0)
  for(j in 1:d){
    lay <- c(lay,rep(j,N[j]))
  }
  
  #allocation and suballocation, each column stores an mcmc output 
  z_post <- matrix(NA,V,n_iter)
  t_post <- matrix(NA,V,n_iter)
  
  #hyperpars
  theta_post <- rep(NA,n_iter)
  theta0_post <- rep(NA,n_iter)
  
  
  # ----------------------------------------------
  # Initialization 
  # ----------------------------------------------
  
  #initialize matrix version of the allocation vector: Z[v,h]=1 if node v is in cluster h, Z[v,h]=0 otherwise
  Z <- vec2mat(z_init)
  
  #initialize list version of suballocation vector: all singletons, coherent with every initialization of the allocations
  T <- vector("list", d)
  for (j in 1:length(N)){
    T[[j]]<-c(1:N[j])
  }
  
  #initialize summaries l and q:
  # matrix l[j,h] = number of subgroups in layer j and cluster h
  # list q[[j]][h,t] = number of nodes in layer j in cluster h and subgroup t (t runs layer-wise)
  l<-matrix(0,d,ncol(Z))
  q <- vector("list", d)
  for (j in 1:d){
    l[j,]  <- colSums(as.matrix(Z[lay==j,]))
    q[[j]] <- t(as.matrix(Z[lay==j,]))
  }  
  
  
  #initialize count matrix of block connections
  temp   <- Y%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  #initialize hyperpars
  theta  <- alpha/lambda
  theta0 <- alpha0/lambda0
  
  
  # ----------------------------------------------
  # Gibbs sampler
  # ----------------------------------------------
  
  for (s in 1:n_iter){
    i<-0
    j<-1
    
    #running across nodes
    for (v in 1:V){
      
      #set layer counter and in-layer node counter
      i<-i+1
      if(i>N[j]){
        i<-1
        j<-j+1}
      
      ## Remove v from Z and l                                                  ##
      ## if the cluster containing node v has no other node, discard it as well ##
      
      #non-empty clusers after removing v
      nonempty_v <- which(colSums(as.matrix(Z[-v,])) > 0)
      if(ncol(Z) > 1){
        
        #reduce Z and l
        Z <- Z[, nonempty_v]
        l <- l[, nonempty_v]
        
        #if it's just one cluster, keep Z and l as matrices
        if (length(nonempty_v)==1){
          Z <- matrix(Z,V,1)
          l <- matrix(l,d,1)}
        
        #reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      
      # Z_v = allocation matrix excluding node v
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y[-v,v])
      
      # h_v = cluster of the removed node, empty if it was a singleton
      h_v         <- which(Z[v,] > 0)
      
      ## Remove node v from q and T                                                ##
      ## if the subgroup containing node v has no other node, discard it as well ##
      
      #matrix version of the suballocation vector of layer j
      Tmat<-vec2mat(T[[j]])
      #non-empty subgroups after removing v
      nonempty_vt <- which(colSums(as.matrix(Tmat[-i,]))>0)
      
      #reduce T and q[[j]]
      if(ncol(Tmat)>1){
        Tmat        <- as.matrix(Tmat[,nonempty_vt])
        q[[j]]      <- as.matrix(q[[j]][nonempty_v,nonempty_vt])
        if(length(nonempty_v)==1){q[[j]] <- t(q[[j]])}
      }
      
      #if the cluster was removed, reduce the rest of q accordingly
      if(length(h_v)!=1){
        for(k in 1:d){
          if(k!=j){
            q[[k]] <- as.matrix(q[[k]][nonempty_v,])
            if(length(nonempty_v)==1)
            {q[[k]] <- t(q[[k]])}
          }
        }
        
      }
      
      ## Compute the count matrices of block connections (and non-connections) after removing v ##
      
      # subgroup of the removed node v
      t_v <- which(Tmat[i,] > 0)
      
      if(length(h_v)==1){# if didn't close a cluster
        resid1       <- matrix(0,H,H)
        resid1[,h_v] <- r_v; resid1[h_v,] <- r_v #residual connections are just r_v
        m            <- m_full - resid1 # compute the m matrix by difference
        if(length(t_v)==1){# if didn't close a subgroup
          q[[j]][h_v,T[[j]][i]] <- q[[j]][h_v,T[[j]][i]]-1 #one less node in the subgroup
        }else{l[j,h_v] <- l[j,h_v]-1} #if closed a subgroup, one less subgroup in the cluster
      }else {m <- m_full} # if closed a cluster (and hence a subgroup), already deleted rows and columns from q, l and m
      
      #rewrite the subgroup vector allocation T from Tmat
      T[[j]] <- Tmat %*% c(1:ncol(Tmat))
      
      #number of subgroups excluding v
      ntab <- ncol(q[[j]])
      
      
      # m_bar = number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      
      #matrix versions of v_minus and r_v
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      ## Compute the probabilities and sample cluster and subgroup for v##
      
      #loglikelihood
      log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      
      #summing predictive probabilities from joint urn scheme
      log_p <- log(joint_HDP(l,q[[j]],N[j],theta,theta0)) + rep(c(log_lhds_old, log_lhd_new),ntab+1)
      p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      
      #sampling allocation and suballocation for v
      new_sample <- matrix(rmultinom(1,1,p),H+1,ntab+1)
      t_new <- which(colSums(new_sample)==1)
      h_new <- which(rowSums(new_sample)==1)
      
      ## Adjust Z, H, r_v, m, T, q and l ##
      
      #write the subllocation in T
      T[[j]][i]<-t_new
      
      #if the old cluster of v still exists, delete the allocation of v in it from Z
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      #if new subgroup
      if(t_new > max(T[[j]][-i])){
        q[[j]] <- cbind(q[[j]],rep(0,H))
        #... and new cluster
        if(h_new == H+1){
          Z               <- cbind(Z,rep(0,V)) #add a column to Z
          Z[v,h_new] <- 1 #allocate v in the new column
          m               <- rbind(cbind(m,0),0) #add a column to m
          r_v             <- crossprod(Z[-v,], Y[-v,v]) #redefine the residual connections
          l               <- cbind(l,rep(0,d)) #add a column to l
          l[j,h_new] <- 1 #count a new subgroup
          for(k in 1:d){
            if(k==j){
              q[[j]] <- rbind(q[[j]],c(rep(0,ntab),1)) #add a column to q[[j]] and count a new node
            }else{q[[k]] <- rbind(q[[k]],rep(0,ncol(q[[k]])))} #add a column to the rest of q
          }
          H <- H + 1 #add a cluster
        }else{#if old cluster
          l[j,h_new]            <- l[j,h_new]+1 #add a new subgroup
          q[[j]][h_new,ntab+1]  <- 1 #count a new node in the new subgroup and old cluster
          Z[v, h_new]           <- 1} #allocate v
      }else{#if old subgroup
        Z[v, h_new]    <- 1  #allocate v
        q[[j]][h_new,t_new] <- q[[j]][h_new,t_new]+1} #add a new node in the old subgroup and old cluster
      
      # update m_full
      resid2              <- matrix(0,H,H)
      resid2[,h_new] <- r_v; resid2[h_new,] <- r_v
      m_full              <- m + resid2
    }
    #end running across nodes
    
    #sample auxiliary variables
    eta0 <- rbeta(1,theta0,sum(l))
    eta  <- sapply(N,function(iter) rbeta(1,1,N))
    
    #sample hyperparameters
    theta0 <- rgamma(1,alpha0+H, rate=(lambda0-log(eta0)))
    theta  <- rgamma(1,alpha+sum(l),rate=(lambda-sum(log(eta))))
    
    #store reordered allocation
    z          <- Z %*% c(1:ncol(Z))
    z_post[,s] <- as.numeric(fct_relabel(fct_inorder(as.factor(z)),~paste0(c(1:length(unique(z))))))
    
    
    #store reordered suballocation
    t <- T
    for(r in 1:d){
      t[[r]] <- as.numeric(fct_relabel(fct_inorder(as.factor(t[[r]])),~paste0(c(1:length(unique(t[[r]]))))))
    }
    t_post[,s] <- list2vec(t,V)
    
    #store hyperpars
    theta_post[s]  <- theta 
    theta0_post[s] <- theta0 
    
    #iteration counter
    if (s%%1000 == 0){print(paste("Iteration:", s))}
  }
  
  # ----------------------------------------------
  # End of Gibbs sampler
  # ----------------------------------------------
  
  #store model specifications and posterior samples in the output list
  pars <- c(a,b,alpha,lambda,alpha0,lambda0)
  model<-list(z_post,t_post,theta_post,theta0_post,pars)
  
  return(model)
}

#EDGE PREDICTION

pred_conn_hyp <- function(Y,lay,newlay,model,nmc){
  #input: Y observed adjacency matrix
  #       lay vector of layer allocations of nodes in Y
  #       newlay vector of layer allocations of new nodes
  #       model output of pex-sbm function
  #       nmc number of MC iterations
  #output: list(adjacency matrix completed with posterior edge probabilities of new nodes, vector of number of new discoveries)
  
  #unpack model characteristics and posterior allocation samples
  N <- dim(Y)[1]
  m <- length(newlay)
  zpost <- model[[1]]
  wpost <- model[[2]]
  thetapost <- model[[3]]
  theta0post <- model[[4]]
  
  ab <- as.numeric(model[[5]][1:2])
  
  nmcmc <- dim(zpost)[[2]]
  
  #check number of mc iterations
  if(nmc>nmcmc){print(paste("Number of MC iterations too large! Using instead number of MCMC iterations:", nmcmc))
    nmc<-nmcmc}
  
  #initialize vector of total and old number of clusters 
  H <- rep(0,nmc)
  Hold <- rep(0,nmc)
  
  #initialize extended adjacency matrix
  Ynew <- cbind(nmc*Y,matrix(0,N,m))
  Ynew <- rbind(Ynew,matrix(0,m,N+m))
  
  #obtain connection probabilities of new nodes for each MC iteration and sum them up
  for (n in 1:nmc){
    
    #consider the (burnin+n)-th joint posterior allocation sample
    z <- zpost[,nmcmc-nmc+n]
    w <- wpost[,nmcmc-nmc+n]
    theta  <- thetapost[nmcmc-nmc+n]
    theta0 <- theta0post[nmcmc-nmc+n]
    laytemp <- lay
    
    #sample sequentially allocations for new nodes via the urn scheme
    for (i in 1:m){
      
      #get summaries for (z_1,...,z_N,...,z_(N+i-1))
      out <- getsummary(z,w,laytemp)
      l <- out[[1]]
      q <- out[[2]]
      ell <- newlay[i] #(N+i)-th node's layer
      Htemp <- max(z) 
      T <- max(w[laytemp==ell]) #number of subgroups in layer ell
      
      #sample jointly allocation and suballocation for (N+i)-th node
      p <- joint_HDP(l,q[[ell]],N,theta,theta0)
      new_sample <- matrix(rmultinom(1,1,p),Htemp+1,T+1)
      wnew <- which(colSums(new_sample)==1)
      znew <- which(rowSums(new_sample)==1)
      
      #update z, w and lay: (N+i)-th node is old now!
      z <- c(z,znew)
      w <- c(w,wnew)
      laytemp <- c(laytemp,newlay[i])
    }
    
    
    #get the matrices M0 and M0bar of counts of edges and non-edges between clusters
    zold <- z[1:N]
    Z    <- vec2mat(zold)
    temp <- Y%*%Z
    M0   <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
    Mtot <- as.matrix(as.numeric(table(zold)))%*%t(as.matrix(as.numeric(table(zold))))
    diag(Mtot) <- as.numeric(table(zold))*(as.numeric(table(zold))-1)/2
    M0bar <- Mtot - M0
    
    H[n] <- max(z)
    Hold[n] <- max(z[1:N])
    M0   <- cbind(M0,matrix(0,Hold[n],H[n]-Hold[n]))
    M0   <- rbind(M0,matrix(0,H[n]-Hold[n],H[n]))
    M0bar   <- cbind(M0bar,matrix(0,Hold[n],H[n]-Hold[n]))
    M0bar   <- rbind(M0bar,matrix(0,H[n]-Hold[n],H[n]))
    
    
    #compute the posterior beta expectation for each extended MC sample 
    for (i in 1:m){
      for(j in 1:(N+i-1)){
        Ynew[N+i,j]<-Ynew[N+i,j]+(ab[1]+M0[z[N+i],z[j]])/(ab[1]+ab[2]+M0[z[N+i],z[j]]+M0bar[z[N+i],z[j]])
        Ynew[j,N+i]<-Ynew[N+i,j]
      }
    }
    if(n%%1000==0){print(paste("Iteration:",n))}
  } 
  return(list(Ynew/nmc,H-Hold))
  
}


####################################################################
################# AUXILIARY FUNCTIONS ##############################
####################################################################

#VECTOR, MATRIX AND LIST TRANSFORMATION

vec2mat <- function(clust_lab){
  # in: vector clust_lab of length V s.t. clust_lab[v]=h if node v is in cluster h
  # out: binary VxH matrix M s.t. M[v,h]=1{node v is in cluster h}
  V <- length(clust_lab)
  H <- max(clust_lab)
  M <- matrix(0,V,H)
  for (v in 1:V){
    M[v,clust_lab[v]] <- 1
  }
  return(M)
}

list2vec <- function(table_list,V){
  # in: list table_list of length d s.t. table_list[[j]][i]=t if the i-th node of layer j is in the t-th subgroup of layer j
  #   : V total number of nodes
  # out: a vector of length V s.t. out[v]=t if node v is in the t-th subgroup of its layer
  d <- length(table_list)
  out <- rep(0,V)
  N <- 0
  cumul <- 0
  for (j in 1:d){
    cumul <- cumul + N
    t <- cumul + 1
    N <- length(table_list[[j]])
    for(i in 1:N){
      out[t] <- table_list[[j]][i]
      t <- t+1
    }
  }
  return(out)
}

vec2list <- function(vec,lay){
  #divides the vector vec in a list of vectors according to the allocation lay
  d <- max(lay)
  out <- vector("list",d)
  cumul<-0
  for (j in 1:d){
    out[[j]]<-vec[lay==j]
  }
  
  return(out)
}

mat2vec <- function(Y){
  #INPUT: a co-clustering matrix (0-1)
  #OUTPUT: the vector of cluster labels
  
  V <- ncol(Y)
  z <- c(1:V)
  for(i in 1:V){
    z[which(Y[i,]==1)]<-z[i]
  }
  return(z)
  
}

#POSTERIOR CO-CLUSTERING MATRIX FROM POSTERIOR SAMPLE

pr_cc <- function(z_post){
  # in: posterior sample of assignments (VxN_iter matrix)
  # out: VxV matrix c with elements c[vu]=fraction of iterations in which v and u are in the same cluster
  V <- nrow(z_post)    
  n_iter <- ncol(z_post)
  c <- matrix(0,V,V)
  for (t in 1:n_iter){
    Z <- vec2mat(z_post[,t])
    c <- c + Z%*%t(Z)
  }
  return(c/n_iter)
}


#SUMMARIES l AND q FROM ALLOCATION AND SUBALLOCATION VECTORS

qfill <- function(Z,W,H){
  #input: Z, W cluster and subgroup allocation vectors in 1 layer, H total number of clusters (not just in the considered layer!)
  #output: vector of the number of nodes in a subgroup AND a cluster, for 1 specific layer
  T <- max(W)
  q <- numeric(0)
  
  for (t in 1:T)
  {q <- cbind(q,replace(numeric(H),unique(Z[W==t]),length(Z[W==t])))}
  return(q)
}

lfill <- function(Z,W,h){
  #input: Z, W cluster and subgroup allocation vectors in 1 layer, h a cluster label
  #output: number of subgroups in h in the layer of interest
  return(length(unique(W[Z==h])))
}

getsummary <- function(Z,W,lay){
  #input: vectors of cluster, subgroup and layer allocations
  #output: list(l,q), l is the matrix l[j,h]=number of subgroups in layer j and cluster h
  #                   q is the list q[[j]][h,t]=number of nodes in layer j, cluster h and subgroup t
  N <- length(Z)
  d <- max(lay)
  
  Wlist <- vec2list(W,lay)
  Zlist <- vec2list(Z,lay)
  
  H     <- max(Z)
  Hlist <- vec2list(rep(H,d),c(1:d))
  
  q <- mapply(qfill,Zlist,Wlist,Hlist, SIMPLIFY=F)
  
  l <- matrix(0,d,H)
  for(h in 1:H){
    l[,h] <- mapply(lfill,Zlist,Wlist,MoreArgs = list(h = h))
  }
  return(list(l,q))
}

############################################################
################# COMPETITORS ##############################
############################################################

####ESBM####

esbm <- function(Y, seed, N_iter, prior, z_init=c(1:nrow(Y)), a=1, b=1,
                 alpha_PY=NA, sigma_PY=NA, beta_DM=NA, H_DM=NA, gamma_GN=NA, x=NULL, alpha_xi=NULL){
  
  # ----------------------------------------------
  # Selection of the prior distribution to be used
  # ----------------------------------------------
  
  if (prior=="DP"){
    urn<-function(v_minus){return(urn_DP(v_minus,alpha_PY))}
  } else if (prior=="PY"){
    urn<-function(v_minus){return(urn_PY(v_minus,alpha_PY,sigma_PY))}
  } else if (prior=="DM"){
    urn<-function(v_minus){return(urn_DM(v_minus,beta_DM,H_DM))}
  } else if (prior=="GN"){
    urn<-function(v_minus){return(urn_GN(v_minus,gamma_GN))}
  } else { 
    stop("Invalid value for prior")  
  }
  
  # ----------------------------------------------
  # Pre-processing of the node attributes
  # ----------------------------------------------
  
  if (!is.null(x)){
    print("Covariates have been provided")
    x <- as.numeric(as.factor(x))
    X <- vec2mat(x)
    if (!is.null(alpha_xi)){
      alpha0 <- sum(alpha_xi)
    } else {
      stop("If covariates x are given, then alpha_xi must be set as well")
    }
  }
  
  # ----------------------------------------------
  # Initialization
  # ----------------------------------------------
  
  set.seed(seed)
  
  V <- nrow(Y)
  
  # cluster assignments are encoded in two equivalent ways:
  # [i] a VxH matrix Z, s.t. Z[v,h]=1{node v in cluster h}, faster to use within each iteration
  Z <- vec2mat(z_init)
  
  # [ii] a vector of length V containing the cluster label for each node, more compact to store;
  # such vectors for all iterations are packed in a VxN_iter matrix z_post, 
  # s.t. z_post[v,t]=h if node v is in cluster h at iteration t
  # Note: matrix z_post takes less memory than a list of N_iter matrices Z
  z_post <- matrix(NA,V,N_iter)
  
  
  # Create the matrix with block connections
  temp   <- Y%*%Z
  m_full <- t(Z)%*%temp - diag(0.5*colSums(temp*Z),ncol(Z))
  
  # ----------------------------------------------
  # Beginning of the Gibbs sampler
  # ----------------------------------------------
  
  for (t in 1:N_iter){
    for (v in 1:V){
      
      # remove empty clusters and
      # if the cluster containing node v has no other node, discard it as well
      if(ncol(Z) > 1){
        nonempty_v <- which(colSums(Z[-v,]) > 0)  
        Z <- Z[, nonempty_v]
        if (length(nonempty_v)==1){Z <- matrix(Z,V,1)}
        
        # Reduce the dimensions of the m_full matrix
        m_full <- matrix(m_full[nonempty_v,nonempty_v],ncol(Z),ncol(Z))
      } 
      
      # H = number of active clusters
      H   <- ncol(Z)
      Z_v <- Z[-v,]
      
      # v_minus = number of nodes in each cluster, excluding node v
      if (H==1){v_minus <- sum(Z[-v])} else {v_minus <- colSums(Z_v)}
      
      # r_v = number of edges from node v to each cluster (no self-loops)
      r_v         <- crossprod(Z_v, Y[-v,v])
      h_v         <- which(Z[v,] > 0)
      
      # Compute the m matrix by difference
      if(length(h_v) == 1){
        resid1       <- matrix(0,H,H)
        resid1[,h_v] <- r_v; resid1[h_v,] <- r_v
        m            <- m_full - resid1
      } else {m <- m_full} # No need to update m in this case
      
      # m_bar = number of non-edges between cluster h and cluster k, excluding node v 
      m_bar   <- matrix(v_minus,H,1)%*%matrix(v_minus,1,H) - diag(0.5*v_minus*(v_minus+1),H) - m
      V_minus <- matrix(1,H,1)%*%matrix(v_minus,1,H)
      R       <- matrix(1,H,1)%*%matrix(r_v,1,H)
      
      # ----------------------------------------------
      # Computing the probabilities
      # ----------------------------------------------
      
      log_lhds_old <- rowSums(lbeta(m + R + a, m_bar + V_minus - R + b) - lbeta(m + a, m_bar + b)) # vector of length H
      log_lhd_new  <- sum(lbeta(r_v + a, v_minus - r_v + b) - lbeta(a, b)) # scalar
      log_addit    <- 0
      
      if(!is.null(x)){
        Vx        <- crossprod(Z_v, X[-v,])
        addit_old <- (Vx[,x[v]] + alpha_xi[x[v]]) / (v_minus+alpha0)
        addit_new <- alpha_xi[x[v]] / alpha0
        log_addit <- log(c(addit_old, addit_new))
      }
      
      # Probabilities
      log_p <- log_addit + log(urn(v_minus)) + c(log_lhds_old, log_lhd_new)
      p     <- exp(log_p - max(log_p)); #p <- p / sum(p)
      
      # ----------------------------------------------
      # Sampling the indicator
      # ----------------------------------------------
      
      new_sample <- which(rmultinom(1,1,p) > 0)
      
      # ----------------------------------------------
      # Adjusting Z, H, r_v and m
      # ----------------------------------------------
      
      if(length(h_v) == 1){Z[v,h_v] <- 0}
      
      # If a new value is sampled, increase the dimension of m_full
      if(new_sample== H+1){
        Z               <- cbind(Z,rep(0,V))
        Z[v,new_sample] <- 1
        m               <- rbind(cbind(m,0),0)
        H               <- H + 1
        r_v             <- crossprod(Z[-v,], Y[-v,v])
      } else {Z[v, new_sample] <- 1}
      
      # Updating m_full
      resid2              <- matrix(0,H,H)
      resid2[,new_sample] <- r_v; resid2[new_sample,] <- r_v
      m_full              <- m + resid2
    }
    
    # store cluster assignments at time t in matrix z_post s.t.
    # z_post[v,t]=h if node v is in cluster h at iteration t
    z_post[,t] <- Z %*% c(1:ncol(Z))
    
    #print(table(z_post[,t])) 
    if (t%%1000 == 0){print(paste("Iteration:", t))}
  }
  return(z_post)
}

urn_DP <- function(v_minus,alpha_PY){
  return(c(v_minus,alpha_PY))
}

urn_PY <- function(v_minus,alpha_PY,sigma_PY){
  H<-length(v_minus)
  return(c(v_minus-sigma_PY,alpha_PY+H*sigma_PY))
}

urn_DM <- function(v_minus,beta_DM,H_DM){
  H<-length(v_minus)
  return(c(v_minus+beta_DM,beta_DM*(H_DM-H)*(H_DM>H)))
}

urn_GN <- function(v_minus,gamma_GN){
  H<-length(v_minus)
  return(c((v_minus+1)*(sum(v_minus)-H+gamma_GN),H^2-H*gamma_GN))
}

sampleLL <- function(memb,Y,a,b){
  # in: vector of cluster labels (memb), VxV adjancency matrix (Y) and hyperparameters beta priors (a,b)
  # out: vector of Bernoulli log-likelihoods for the edges under ESBM (conditioned on memb and on block-probabilities)
  
  z <- vec2mat(memb)
  H <- ncol(z)
  V <- dim(Y)[1]
  
  M <- t(z)%*%Y%*%z
  diag(M) <- diag(M)/2
  Tot <- t(z)%*%matrix(1,V,V)%*%z
  diag(Tot) <- (diag(Tot)-table(memb))/2
  Mbar <- Tot-M
  a_n <- lowerTriangle(M,diag=TRUE)+a
  b_bar_n <- lowerTriangle(Mbar,diag=TRUE)+b
  
  theta <- rbeta(length(a_n),a_n,b_bar_n)
  Theta <- matrix(0,H,H)
  Theta[lower.tri(Theta,diag=TRUE)] <- theta
  Theta <- Theta+t(Theta)
  diag(Theta) <- diag(Theta)/2
  edge_prob <- z%*%Theta%*%t(z)
  
  LL <- dbinom(lowerTriangle(Y,diag=FALSE), size=1, prob=lowerTriangle(edge_prob,diag=FALSE),log=TRUE)
  
  return(LL)
}