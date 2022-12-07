\name{StatComp22083-package}
\alias{StatComp22083-package}
\alias{StatComp22083}
\docType{package}
\title{
  To compute item response probabilities and estimating parameters based on 2PL IRT model.
}
\description{
  This package is created for computing item response probabilities and estimating parameters based on 2PL(discrimination and difficulty parameters) IRT model. Through MCMC method(using Gibbs update), we can obtain estimates of these two parameters and ability parameter of students.
}
\details{
   Use 'estimateR' function to compute item response probabilities and estimating parameters based on 2PL IRT model.
}
\author{
Litian Xu

Maintainer: Litian Xu <xlt@mail.ustc.edu.cn>
}
\keyword{ package }
\examples{
  \dontrun{
    library(ltm)
    estimateR(t(LSAT))   
  }
}
