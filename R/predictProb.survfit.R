predictProb.survfit <- function (object, newdata, times, train.data=NULL, ...){
  km.pred <- matrix(summary(object, times=times)$surv, nrow(newdata), length(times),byrow=TRUE)  
}