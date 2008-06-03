\name{extract.fun}
\alias{extract.fun}
\title{Extract functions and libraries to be loaded onto a compute cluster}
\description{
Automatic extraction of functions and libraries employed passed functions. Designed for \code{peperr} call, see Details section there.
}
\usage{
extract.fun(funs = NULL)
}
\arguments{
\item{funs}{list of function names.}
}
\value{
list containing
  \item{packages}{vector containing quoted names of libraries}
  \item{functions}{vector containing quoted names of functions}

}
\details{
This function is necessary for compute cluster situations where for computation on nodes required functions and libraries have to be loaded explicitly on each node. Avoids loading of whole global environment which might include the unnecessary loading of huge data sets. 
}
\seealso{\code{peperr}}
\examples{
# 1. Simplified example for illustration
\dontrun{
library(CoxBoost)
a <- function(z){
# some calculation
}

b<-function(){
# some other calculation
x <- cv.CoxBoost(values)
y <- a(z)
}

# list with packages and functions required for b:
extract.fun(list(b))

# 2. As called by default in peperr example
extract.fun(list(fit.CoxBoost, aggregation.pmpec))
}
}
\keyword{models} \keyword{regression} \keyword{survival}