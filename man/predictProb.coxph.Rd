\name{predictProb.coxph}
\alias{predictProb.coxph}
\title{Extract predicted survival probabilities from a coxph object}
\description{
Extracts predicted survival probabilities for survival models fitted by Cox proportional hazards model, providing an interface as required by \code{pmpec}. 
}
\usage{
predictProb.coxph(object, newdata, times, train.data,  ...)
}
\arguments{
\item{object}{a fitted model of class \code{coxph}.}
\item{newdata}{data frame containing observations for which predictions are desired.}
\item{times}{vector of evaluation time points.}
\item{train.data}{not used.}
\item{\dots}{additional arguments, currently not used.}
}
\value{
Matrix with probabilities for each evaluation time point in \code{times}(columns) and each new observation (rows). 
}
\keyword{models} \keyword{regression} \keyword{survival}