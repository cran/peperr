\name{predictProb.CoxBoost}
\alias{predictProb.CoxBoost}
\title{Extract predicted survival probabilities from a CoxBoost fit}
\description{
Extracts predicted survival probabilities from survival model fitted by CoxBoost, providing an interface as required by \code{pmpec}. 
}
\usage{
predictProb.CoxBoost(object, newdata, times, train.data=NULL, complexity, ...)
}
\arguments{
\item{object}{a fitted model of class \code{CoxBoost}.}
\item{newdata}{data frame containing observations for which predictions are desired.}
\item{times}{vector of evaluation time points.}
\item{train.data}{not used.}
\item{complexity}{complexity value.}
\item{\dots}{additional arguments, currently not used.}
}

\value{
Matrix with probabilities for each evaluation time point in \code{times} (columns) and each new observation (rows). 
}
\keyword{models} \keyword{regression} \keyword{survival}