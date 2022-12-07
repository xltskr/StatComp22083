#' @title Update theta
#' @description We use gibbs sampler to obtain estimates of parameters, refer to the book Item Response Theory Parameter Estimation Techniques, 2rd edition. During samples generation, we need to update theta,a,b respectively. The following three functions are displayed.Function written for updating ability parameter 'theta'.
#' @param a discriminationparameter of an item
#' @param b difficulty parameters of an item
#' @param theta ability parameter of a participant/subject/student
#' @param x subject's response to an item
#' @param sigma standard error of random normal distribution that is set in advance
#' @export
update_theta <- function(theta,a,b,x,sigma){
  theta_candidate <- rnorm(length(theta),theta,sigma)
  rho <- exp(apply(loglik(theta_candidate,a,b,x),2,sum)
             - apply(loglik(theta,a,b,x),2,sum) + dnorm(theta_candidate,0,1,log=TRUE)
             - dnorm(theta,0,1,log=TRUE))
  ind_update <- 1*(rho > runif(length(theta),0,1))
  theta <- ind_update*theta_candidate + (1-ind_update)*theta
  return (theta)
}