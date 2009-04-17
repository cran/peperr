\name{fit.rsf_mtry}
\alias{fit.rsf_mtry}
\title{Interface function for fitting a random survival forest}
\description{
Interface for fitting a random survival forest by \code{rsf}, conforming to the requirements for argument \code{fit.fun} in \code{peperr} call. 
}
\usage{
fit.rsf_mtry(response, x, cplx, ...)
}
\arguments{
  \item{response}{a survival object (with \code{Surv(time, status)}).}
  \item{x}{\code{n*p} matrix of covariates.}
  \item{cplx}{Value of mtry, the number of randomly selected variables at each split.}
  \item{\dots}{additional arguments passed to \code{rsf} call.}
}
\value{
rsf object
}
\details{
Function is basically a wrapper around \code{rsf} of package \pkg{randomSurvivalForest}. Variable mtry, i.e.\ the number of variables selected randomly at each split, is taken as complexity value.
}
\seealso{ \code{peperr}, \code{\link[randomSurvivalForest]{rsf}}}
\keyword{models} \keyword{regression} \keyword{survival}