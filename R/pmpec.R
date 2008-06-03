pmpec <- function(object, data, times, model.args=NULL, type=c("PErr","NoInf"), 
                  external.time=NULL, external.status=NULL) {
    require(survival)

    type <- match.arg(type)

    fit.time <- data$time
    fit.status <- data$status
    
    if (!is.null(external.time) && !is.null(external.status)) {
        fit.time <- external.time
        fit.status <- external.status
    }

    censdata <- data.frame(time=fit.time,status=ifelse(fit.status == 0,1,0))
    cens.fit <- survival::survfit(Surv(time,status) ~ 1,data=censdata)

    eval.cens.prob <- summary(cens.fit,times=times)$surv
    
    if (is.null(external.time) || is.null(external.status)) {
        sortsurv <- summary(cens.fit,times=sort(data$time))$surv
        invcensprob <- rep(NA,length(data$status))
        invcensprob[order(data$time)] <- c(1,sortsurv[1:length(sortsurv)-1])
    } else {
        sort.time <- sort(fit.time)
        sortsurv <- c(1,summary(cens.fit,times=sort.time)$surv)
        sort.time <- c(-999,sort.time)
        sortsurv <- c(1,sortsurv)
        invcensprob <- unlist(lapply(data$time,function(actual.time) sortsurv[which(sort.time >= actual.time)[1]-1]))
    }

    mod.time <- data$time
    if (any(data$status > 1)) {
        mod.time[data$status > 1] <- max(data$time,times) + 1
    }

    status.mat <- matrix(times,length(data$time),length(times),byrow=TRUE) < mod.time

    probmat <- do.call("predictProb", c(list(object=object,newdata=data,times=times,learn.data=data),model.args))

    weightmat <- t(t(status.mat) / eval.cens.prob) + (1 - status.mat) * matrix(data$status,length(data$status),length(times)) * matrix(1/invcensprob,length(data$status),length(times))

    if (type == "PErr") {
        return(apply(weightmat*(status.mat - probmat)^2,2,mean))
    } else {
        
        #   vectorized version (that actually seems to be slower!)
        # return(apply(matrix(as.vector(t(weightmat)[,rep(1:nrow(data),rep(nrow(data),nrow(data)))])*(
        #                     as.vector(t(status.mat)[,rep(1:nrow(data),rep(nrow(data),nrow(data)))]) - as.vector(t(probmat)))^2
        #                     ,ncol=length(times),byrow=TRUE),2,mean))
        
        res <- .C(noinf,
                  as.double(as.vector(weightmat)),
                  as.integer(as.vector(status.mat)),
                  as.double(as.vector(probmat)),
                  as.integer(nrow(data)),
                  as.integer(length(times)),
                  err=double(length(times))
                  )
        
        return(res$err)
                
        # t.status.mat  <- t(status.mat)
        # t.weightmat <- t(weightmat)        
        # errmat <- matrix(unlist(lapply(1:nrow(probmat),function(i) apply(t.weightmat*(t.status.mat - probmat[i,])^2,1,mean))),nrow(probmat),length(times),byrow=TRUE)
        # return(apply(errmat,2,mean))
    }
}

