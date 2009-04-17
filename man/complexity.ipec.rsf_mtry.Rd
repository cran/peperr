\name{complexity.ipec.rsf_mtry}
\alias{complexity.ipec.rsf_mtry}
\title{Interface function for complexity selection for random survival forest via integrated prediction error curve and the bootstrap}
\description{
Determines the optimal value for the number of candidate values at each split for a survival model fitted by random survival forest via integrated prediction error curve (IPEC) estimates, conforming to the calling convention required by argument \code{complexity} in \code{peperr} call. 
}
\usage{
complexity.ipec.rsf_mtry(response, x, boot.n.c = 10, mtry, 
   eval.times = NULL, full.data, relative=TRUE, ...)
}
\arguments{
  \item{response}{a survival object (with \code{Surv(time, status)}).}
  \item{x}{\code{n*p} matrix of covariates.}
  \item{boot.n.c}{number of bootstrap samples.}
  \item{mtry}{vector of potential values for parameter mtry, i.e. number of randomly selected variables at each split.}
  \item{eval.times}{vector of evaluation time points.}
  \item{full.data}{Data frame containing response and covariates of the full data set.}
  \item{relative}{Should relative IPEC be used for evaluation?} 
  \item{\dots}{additional arguments passed to \code{rsf} call.}
}
\value{
Scalar value giving the optimal value for parameter mtry.
}
\details{
Plotting the .632+ estimator for each time point given in \code{eval.times} results in a prediction error curve. A summary measure can be obtained by integrating over time. \code{complexity.ipec.rsf_mtry} uses a Lebesgue-like integral taking Kaplan-Meier estimates as weights. The optimal value for parameter mtry out of the set of passed values, for which the minimal IPEC is obtained, is returned.
}
\seealso{ \code{peperr}, \code{\link[randomSurvivalForest]{rsf}} }
\keyword{models} \keyword{regression} \keyword{survival}
