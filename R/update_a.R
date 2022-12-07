#' @title Update a
#' @description Function written for updating discrimination parameter 'a'.
#' @param a discriminationparameter of an item
#' @param b difficulty parameters of an item
#' @param theta ability parameter of a participant/subject/student
#' @param x subject's response to an item
#' @param sigma standard error of random normal distribution that is set in advance
#' @export
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