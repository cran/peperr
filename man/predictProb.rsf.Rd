\name{predictProb.rsf}
\alias{predictProb.rsf}
\title{Extract predicted survival probabilities from a random survival forest fit}
\description{
Extracts predicted survival probabilities from survival model fitted by randomSurvivalForest function \code{rsf}, providing an interface as required by \code{pmpec}. 
}
\usage{
predictProb.rsf(object, response, x, times, ...)
}
\arguments{
\item{object}{a fitted model of class \code{rsf}.}
\item{response}{survival object (with \code{Surv(time, status)}, where time is an \code{n}-vector of censored survival times and status an \code{n}-vector containing survival status, coded with 0 and 1.}
\item{x}{\code{n*p} matrix of covariates.}
\item{times}{vector of evaluation time points.}
\item{\dots}{additional arguments, currently not used.}
}

\value{
Matrix with probabilities for each evaluation time point in \code{times} (columns) and each new observation (rows). 
}
\keyword{models} \keyword{regression} \keyword{survival}