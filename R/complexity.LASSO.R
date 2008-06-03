complexity.LASSO <- function(response, x, full.data, ...){
   require(penalized)
   lambda <- optL1(response=response, penalized=x, ...)$lambda
   lambda
}