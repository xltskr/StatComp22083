#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::export]]
NumericMatrix gibbsC(double mu1,double mu2,double sigma1,double sigma2,double rho,int N) {
  NumericMatrix mat(N, 2);
  double s1=sqrt(1-pow(rho,2))*sigma1;
  double s2=sqrt(1-pow(rho,2))*sigma2;
  mat(0,0)=mu1;
  mat(0,1)=mu2;
  for(int i = 1; i < N; i++) {
    double x2=mat(i-1,1);
    double m1=mu1+rho*(x2-mu2)*sigma1/sigma2;
    mat(i,0)=rnorm(1,m1,s1)[0];
    double x1=mat(i,0);
    double m2=mu2+rho*(x1-mu1)*sigma2/sigma1;
    mat(i,1)=rnorm(1,m2,s2)[0];
  }
  return(mat);
}