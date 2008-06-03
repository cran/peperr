predictProb.CoxBoost <- function(object, newdata, times, train.data=NULL, complexity, ...){
  predict(object, type="risk", newdata=as.matrix(newdata[,!(names(newdata) %in% c("time","status")),drop=FALSE]),
      times=times, at.step=complexity)
}