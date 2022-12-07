#' @title Compute the log-likelihood function of item response probability.
#' @description Due to distichous x, f(x|p)=p^x(1-p)^{1-x}. The log-likelihood function is x*log(p)+(1-x)log(1-p).
#' @param a discriminationparameter of an item
#' @param b difficulty parameters of an item
#' @param theta ability parameter of a participant/subject/student
#' @param x subject's response to an item 
#' @examples
#' \dontrun{
#' loglik(2,0.5,1,1)
#' }
#' @export
loglik <- function(theta,a,b,x){
  p <- irt_2pl(a,b,theta)
  logp <- x*log(p) + (1-x)*log(1-p)
  return (logp)
}