#' @title Parameter estimates of 2PL IRT model by MCMC method based on Gibbs sampler.
#' @description A Gibbs sampler using R. Due to existence of burn-in period, we take part of sample values and take the average. What's more, the following function permits several chains to ensure convergence.
#' @param x dataset include dichotomous responses to items of subjects 
#' @param n_chain number of chain
#' @param burn burnin length
#' @param n_itr length of each chain
#' @return estimates of theta, a and b
#' @import Rcpp ltm MASS bootstrap boot microbenchmark DAAG scales latex2exp knitr pander
#' @importFrom stats rnorm dlnorm runif dnorm rlnorm
#' @examples
#' \dontrun{
#' library(ltm)
#' est<-estimateR(t(LSAT))
#' est$theta
#' est$a
#' est$b
#' }
#' @export
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