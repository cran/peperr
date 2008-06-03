\name{predictProb}
\alias{predictProb}
\title{Generic function for extracting predicted survival probabilities}
\description{
Generic function for extraction of predicted survival probabilities from a fitted survival model conforming to the interface required by \code{pmpec}.
}
\usage{
predictProb(object, ...)
}
\arguments{
\item{object}{a fitted survival model.}
\item{\dots}{additional arguments.}
}
\details{
\code{pmpec} requires a \code{predictProb.class} function for each model fit of class \code{class}. It extracts the predicted probability of survival from this model.

See existing \code{predictProb} functions, at the time \code{predictProb.CoxBoost} and \code{predictProb.survfit}, or \code{predictProb} functions provided by package \pkg{pec}. 
}

\value{
Matrix with predicted probabilities for each evaluation time point in \code{times} (columns) and each new observation (rows). 
}
\keyword{models} \keyword{regression} \keyword{survival}