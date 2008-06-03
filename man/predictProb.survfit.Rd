\name{predictProb.survfit}
\alias{predictProb.survfit}
\title{Extract predicted survival probabilities from a survfit object}
\description{
Extracts predicted survival probabilities for survival models fitted by \code{survfit}, providing an interface as required by \code{pmpec}. 
}
\usage{
predictProb.survfit(object, newdata, times, train.data,  ...)
}
\arguments{
\item{object}{a fitted model of class \code{survfit}.}
\item{newdata}{data frame containing observations for which predictions are desired.}
\item{times}{vector of evaluation time points.}
\item{train.data}{not used.}
\item{\dots}{additional arguments, currently not used.}
}
\value{
Matrix with probabilities for each evaluation time point in \code{times}(columns) and each new observation (rows). 
}
\keyword{models} \keyword{regression} \keyword{survival}