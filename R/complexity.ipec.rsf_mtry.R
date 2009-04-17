complexity.ipec.rsf_mtry <-
function(response, x, boot.n.c=10, mtry, 
   eval.times=NULL, full.data, relative=TRUE,  ...){
   require(randomSurvivalForest)
   require(survival)
   actual.data.c <- as.data.frame(x)
   xnames <- names(actual.data.c)
   time <- response[,"time"]
   status <- response[,"status"]
   actual.data.c$time <- time
   actual.data.c$status <- status

   n.c <- length(time)
   boot.index.c <- matrix(sapply(1:boot.n.c,function(b){sample(1:n.c,replace=TRUE)}), 
      byrow=TRUE, nrow=boot.n.c, ncol=n.c)
   not.in.sample.c <- list()
   for (i in 1:boot.n.c){ 
      not.in.sample.c[[i]] <- (1:n.c)[-unique(boot.index.c[i,])]
   }

   uncens <- which(status == 1)
   if (is.null(eval.times)){
      if(length(unique(time[uncens]))<100){
         w.eval.times.c <- c(0,sort(time[uncens]))
      } else {
      w.quantile <- 90
      space <- round(length(time[uncens])/100)
      index <- (1:w.quantile)*space
      w.eval.times.c <- c(0,sort(time[uncens])[index[index<length(time[uncens])]])
      }
   } else {
      w.eval.times.c <- sort(eval.times)
   }

   w.eval.times.c <- unique(w.eval.times.c)

   km.fit <- survfit(Surv(time, status)~1, data=actual.data.c)
   km.pred <- summary(object=km.fit, times=w.eval.times.c)$surv
   km.weight <- -1*diff(km.pred)

   fullforest <- list()
   counter <- 1
   for (i in mtry){
   fullforest[[counter]] <- rsf(Survrsf(time, status)~., data=actual.data.c, mtry=i, ...)
   counter <- counter + 1
   }

   w.full.apparent <- matrix(NA,nrow=length(mtry), ncol=length(w.eval.times.c))
   w.noinf.error <- matrix(NA,nrow=length(mtry), ncol=length(w.eval.times.c))

   for (m in 1:length(mtry)){
      pe.w.full.apparent <- pmpec(object=fullforest[[m]],
         data=actual.data.c, times=w.eval.times.c,
         #model.args=list(complexity=m), überflüssig?
         external.time=full.data$time, external.status=full.data$status)
      w.full.apparent[m,] <- pe.w.full.apparent

      pe.w.noinf.error <- pmpec(object=fullforest[[m]],
         data=actual.data.c, times=w.eval.times.c,
         #model.args=list(complexity=m), überflüssig
         external.time=full.data$time, external.status=full.data$status, type="NoInf")
      w.noinf.error[m,] <- pe.w.noinf.error
   }

   w.boot.error.wo <- array(dim=c(length(mtry), boot.n.c, length(w.eval.times.c)))

   for (actual.boot in 1:boot.n.c) {
      boot.fit <- list()
      boot.data <- as.data.frame(x[boot.index.c[actual.boot,],])
      boot.data$time <- time[boot.index.c[actual.boot,]]
      boot.data$status <- status[boot.index.c[actual.boot,]]
      counter <- 1
      for (i in mtry){ 
         boot.fit[[counter]] <- rsf(Survrsf(time, status)~., data=boot.data, mtry=i, ...)
      counter <- counter + 1
      }
     
      for (m in 1:length(mtry)){
         w.pec.boot <- pmpec(object=boot.fit[[m]],
            data=actual.data.c[not.in.sample.c[[actual.boot]],], 
            times=w.eval.times.c, 
        #    model.args=list(complexity=m), überflüssig?
            external.time=full.data$time, external.status=full.data$status)
         w.boot.error.wo[m,actual.boot,] <- w.pec.boot
      }
   }

   if(relative){
      for (i in 1:boot.n.c) w.boot.error.wo[,i,] <- t(t(w.boot.error.wo[,i,]) - w.boot.error.wo[1,i,])
         w.noinf.error <- t(t(w.noinf.error) - w.noinf.error[1,])
         w.full.apparent <- t(t(w.full.apparent) - w.full.apparent[1,])
   }

   w.mean.boot.error.wo <- apply(w.boot.error.wo,c(1,3),mean,na.rm=TRUE)
   w.boot632p.error.wo <- matrix(NA, nrow=length(mtry), ncol=length(w.eval.times.c))

   for (m in 1:length(mtry)) {

      w.relative.overfit.wo <- ifelse(w.noinf.error[m,] > w.full.apparent[m,],
         (ifelse(w.mean.boot.error.wo[m,] < w.noinf.error[m,], 
         w.mean.boot.error.wo[m,], w.noinf.error[m,]) - w.full.apparent[m,])/
         (w.noinf.error[m,] - w.full.apparent[m,]), 0)
      w.weights <- .632/(1-.368*w.relative.overfit.wo)

      w.boot632p.error.wo[m,] <- (1-w.weights)*w.full.apparent[m,] + 
         w.weights*ifelse(w.mean.boot.error.wo[m,] < w.noinf.error[m,], 
         w.mean.boot.error.wo[m,], w.noinf.error[m,])
   }

   Lint.boot632p.error <- apply(t(w.boot632p.error.wo[,1:(length(km.weight))])*km.weight,2,sum)
  
   mtry[which.min(Lint.boot632p.error)]
}
