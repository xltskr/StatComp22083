#include <Rcpp.h>
using namespace Rcpp;

//' @title Compute the log-likelihood function  using Rcpp
//' @param theta ability parameter of a participant/subject/student
//' @param a discriminationparameter of an item
//' @param b difficulty parameters of an item
//' @param x subject's response to an item 
//' @examples
//' \dontrun{
//' lp<-loglikC(2,0.5,1,1)
//' }
//' @useDynLib StatComp22083
//' @export
// [[Rcpp::export]]
double loglikC(double theta, double a, double b, int x){
  double p = 1 / (1 + exp(a*(b-theta)));
  double logp = x*log(p) + (1-x)*log(1-p);
  return logp;
}