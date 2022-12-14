% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglik.R
\name{loglik}
\alias{loglik}
\title{Compute the log-likelihood function of item response probability.}
\usage{
loglik(theta, a, b, x)
}
\arguments{
\item{theta}{ability parameter of a participant/subject/student}

\item{a}{discriminationparameter of an item}

\item{b}{difficulty parameters of an item}

\item{x}{subject's response to an item}
}
\description{
Due to distichous x, f(x|p)=p^x(1-p)^{1-x}. The log-likelihood function is x*log(p)+(1-x)log(1-p).
}
\examples{
\dontrun{
loglik(2,0.5,1,1)
}
}
