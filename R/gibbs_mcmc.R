#' @title Posterior estimate
#' @description Function to obtain posterior estimates of theta,a and b.
#' @param x dataset include dichotomous responses to items of subjects 
#' @param n_itr length of each chain
#' @return posterior estimates of theta,a and b for n_itr iterations.
#' @examples
#' \dontrun{
#' librar(ltm)
#' gibbs_mcmc(t(LSAT),3e3)
#' }
#' @export
gibbs_mcmc <- function(x,n_itr){
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