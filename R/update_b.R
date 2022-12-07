#' @title Update b
#' @description Function written for updating difficulty parameter 'b'.
#' @param a discriminationparameter of an item
#' @param b difficulty parameters of an item
#' @param theta ability parameter of a participant/subject/student
#' @param x subject's response to an item
#' @param sigma standard error of random normal distribution that is set in advance
#' @export
update_b <- function(b,theta,a,x,sigma){
  b_candidate <- rnorm(length(b),b,sigma)
  rho <- exp(apply(loglik(theta,a,b_candidate,x),1,sum)
             - apply(loglik(theta,a,b,x),1,sum) + dnorm(b_candidate,0,1,log=TRUE)
             - dnorm(b,0,1,log=TRUE))
  ind_update <- 1*(rho > runif(length(b),0,1))
  b <- ind_update*b_candidate + (1-ind_update)*b
  return (b)
}