#' @title compute item response probability
#' @description Under 2PL IRT model, we consider responses to item are distichous, thus, the probability that one subject with ability theta response correctly to an item with dicrimination 'a' and difficulty 'b' is P(X=1)=expit(a(theta-b))
#' @param a discrimination parameter of an item
#' @param b difficulty parameters of an item
#' @param theta ability parameter of a participant/subject/student
#' @examples
#' \dontrun{
#' irt_2pl(1,0.5,2)
#' }
#' @export
irt_2pl <- function(a,b,theta){
  #compute item response probabilities based on 2PL. The parameterization
  # of the item response function follows the classic 2PL (discrimination and
  # difficulty parameters). The function is vectorized.
  p <- 1 / (1 + exp(a*outer(b,theta,"-")))
  return (p)
}