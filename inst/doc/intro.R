## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----eval=FALSE---------------------------------------------------------------
#  loglikC(2,0.5,1,1)

## -----------------------------------------------------------------------------
irt_2pl <- function(a,b,theta){
  p <- 1 / (1 + exp(a*outer(b,theta,"-")))
  return (p)
}

## -----------------------------------------------------------------------------
loglik <- function(theta,a,b,x){
  p <- irt_2pl(a,b,theta)
  logp <- x*log(p) + (1-x)*log(1-p)
  return (logp)
}

## -----------------------------------------------------------------------------
update_theta <- function(theta,a,b,x,sigma){
  theta_candidate <- rnorm(length(theta),theta,sigma)
  rho <- exp(apply(loglik(theta_candidate,a,b,x),2,sum)
  - apply(loglik(theta,a,b,x),2,sum) + dnorm(theta_candidate,0,1,log=TRUE)
  - dnorm(theta,0,1,log=TRUE))
  ind_update <- 1*(rho > runif(length(theta),0,1))
  theta <- ind_update*theta_candidate + (1-ind_update)*theta
  return (theta)
}

## -----------------------------------------------------------------------------
update_a <- function(a,theta,b,x,sigma){
  a_candidate <- rnorm(length(a),a,sigma)
  if(any(a_candidate < 0)){
    a_candidate[a_candidate < 0] <- a[a_candidate < 0] #automatically reject negative candidates
  }
  rho <- exp(apply(loglik(theta,a_candidate,b,x),1,sum)
    - apply(loglik(theta,a,b,x),1,sum) + dlnorm(a_candidate,0,0.5,log=TRUE)
    - dlnorm(a,0,0.5,log=TRUE))
  ind_update <- 1*(rho > runif(length(a),0,1))
  a <- ind_update*a_candidate + (1-ind_update)*a
  return (a)
}

## -----------------------------------------------------------------------------
update_b <- function(b,theta,a,x,sigma){
  b_candidate <- rnorm(length(b),b,sigma)
  rho <- exp(apply(loglik(theta,a,b_candidate,x),1,sum)
    - apply(loglik(theta,a,b,x),1,sum) + dnorm(b_candidate,0,1,log=TRUE)
    - dnorm(b,0,1,log=TRUE))
  ind_update <- 1*(rho > runif(length(b),0,1))
  b <- ind_update*b_candidate + (1-ind_update)*b
  return (b)
}

## -----------------------------------------------------------------------------
gibbs_mcmc <- function(x,n_itr=3e3){
  N <- dim(x)[2]
  J <- dim(x)[1]
  theta_post_mat <- matrix(NA,nrow=n_itr,ncol=N)
  a_post_mat <- matrix(NA,nrow=n_itr,ncol=J)
  b_post_mat <- matrix(NA,nrow=n_itr,ncol=J)
  theta_post_mat[1,] <- rnorm(N,0,1)
  a_post_mat[1,] <- rlnorm(J,0,0.5)
  b_post_mat[1,] <- rnorm(J,0,1)
  for (i in 2:n_itr){
    theta_post_mat[i,] <- update_theta(theta_post_mat[i-1,],a_post_mat[i-1,],b_post_mat[i-1,],x,0.1)
    a_post_mat[i,] <- update_a(a_post_mat[i-1,],theta_post_mat[i,],b_post_mat[i-1,],x,0.1)
    b_post_mat[i,] <- update_b(b_post_mat[i-1,],theta_post_mat[i,],a_post_mat[i,],x,0.1)
  }
  out <- list(theta_post_mat,a_post_mat,b_post_mat)
  return (out)
}

## -----------------------------------------------------------------------------
estimateR<-function(x,n_chain=3,burn=1e3,n_itr=3e3){
  d<-dim(x)
  ns<-n_itr-burn
  theta_post<-matrix(nrow=ns*n_chain,ncol=d[2])
  a_post<-matrix(nrow=ns*n_chain,ncol=d[1])
  b_post<-matrix(nrow=ns*n_chain,ncol=d[1])
  for(i in 1:n_chain){
    chain<-gibbs_mcmc(x,n_itr)
    theta_post[((ns*(i-1)+1):(ns*i)),]<-chain[[1]][(burn+1):n_itr,]
    a_post[((ns*(i-1)+1):(ns*i)),]<-chain[[2]][(burn+1):n_itr,]
    b_post[((ns*(i-1)+1):(ns*i)),]<-chain[[3]][(burn+1):n_itr,]
  }
  theta_eap <- apply(theta_post,2,mean)
  a_eap <- apply(a_post,2,mean)
  b_eap <- apply(b_post,2,mean)
  return(list(theta=theta_eap,a=a_eap,b=b_eap))
}

## ----eval=FALSE---------------------------------------------------------------
#  library(ltm)
#  # take the dataset 'LSAT' in package 'ltm' for example
#  est<-estimateR(t(LSAT))
#  est$theta
#  est$a
#  est$b

