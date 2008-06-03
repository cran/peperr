`peperr` <-
function(response, x, 
   indices=NULL,
   fit.fun, complexity=NULL, args.fit=NULL, args.complexity=NULL, 
   parallel=FALSE, cpus=2, noclusterstart=FALSE, noclusterstop=FALSE,
   aggregation.fun=NULL, args.aggregation=NULL,
   load.list=extract.fun(list(fit.fun, complexity, aggregation.fun)), 
   load.vars=NULL, load.all=FALSE,
   trace=FALSE, debug=FALSE, peperr.lib.loc=NULL)
{
   require(survival)
   binary <- FALSE
   if(is.null(aggregation.fun)){
      if(is.Surv(response)){
         aggregation.fun <- aggregation.pmpec
      } else {
         if (is.vector(response)&& length(unique(response))<3){
            aggregation.fun <- aggregation.brier
         } else {
            stop("Please specify argument 'aggregation.fun' according to 'response'")  
         }
      }
   }

   if(is.null(indices)){
      if(nrow(x)<ncol(x)){
         indices <- resample.indices(n=nrow(x), method="sub632", sample.n=500)
      } else {
         indices <- resample.indices(n=nrow(x), method="boot", sample.n=500)
      }
   }

   if(parallel==TRUE){
      sfInit(parallel=TRUE, cpus=cpus, nostart=noclusterstart)
      } else { 
      sfInit(nostart=noclusterstart)
   }
   sfLibrary(peperr, lib.loc=peperr.lib.loc)

   if(load.all==TRUE){
      sfExportAll()
   } else {
      if (!is.null(load.list)){
         for (i in seq(along=load.list$packages)){
            try(eval(call("sfLibrary", load.list$packages[i], 
               character.only=TRUE)), silent=!debug)
         }
         for (i in seq(along=load.list$functions)){

            try(eval(call("sfExport", load.list$functions[i], debug=debug)), silent=!debug)
         }
      }
      if (!is.null(load.vars)){
         sfExport(names(load.vars))
      }
   }

   actual.data <- as.data.frame(x)
   if (is.Surv(response)){
      #xnames <- names(actual.data) 
      time <- response[,"time"]
      status <- response[, "status"]
      actual.data$time <- time
      actual.data$status <- status
      km.fit <- survfit(Surv(time, status)~1, data=actual.data)
      null.model <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
         response=response, x=x, model=km.fit), args.aggregation))
   } else {
      if (is.vector(response)&& length(unique(response))<3){
         binary <- TRUE
         actual.data$response <- response
         logreg.fit <- glm(formula=response~1, data=actual.data, family=binomial())
         null.model <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
            response=response, x=x, model=logreg.fit), args.aggregation))
      }
   }
   dimnames(null.model) <- NULL
   attr(null.model, "addattr") <- attributes(null.model)$addattr
   fullsample.attr <- attr(null.model, "addattr")

   fullsample.index <- resample.indices(n=nrow(x), method="no")
   sample.index <- indices$sample.index
   not.in.sample <- indices$not.in.sample
   sample.index.full <- indices$sample.index
   not.in.sample.full <- indices$not.in.sample
   sample.n <- length(sample.index)+1
   sample.index.full[[sample.n]] <- fullsample.index$sample.index[[1]]
   not.in.sample.full[[sample.n]] <- fullsample.index$not.in.sample[[1]]

   if (is.function(complexity)){
      cplx <- do.call("complexity", c(list(response=response, x=x, full.data=actual.data), 
         args.complexity))
   } else {
      if (is.vector(complexity)){
         cplx <- complexity
      } else cplx <- 0
   }

#    full.apparent <- c()
    noinf.error <- c()
# 
   if (is.list(cplx)){
      for (i in 1:length(cplx[[1]])){
          list.cplx <- lapply(cplx, function(arg) arg[i])
          fullmodel <- do.call("fit.fun", c(list(response=response, x=x, 
            cplx=list.cplx), args.fit))
     
#          full.apparent.i <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
#             response=response, x=x, model=fullmodel, cplx=list.cplx), args.aggregation))	
#          full.apparent <- rbind(full.apparent, full.apparent.i)
         noinf.error.i <- do.call("aggregation.fun", c(list(full.data=actual.data, type="noinf", 
               response=response , x=x, model=fullmodel, cplx=list.cplx), args.aggregation))	
         noinf.error <- rbind(noinf.error, noinf.error.i)	
      }
   } else {
      if (is.vector(cplx)){
         for (i in 1:length(cplx)){
            fullmodel <- do.call("fit.fun", c(list(response=response, x=x, 
              cplx=cplx[i]), args.fit))
    
#             full.apparent.i <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
#                response=response, x=x, model=fullmodel, cplx=cplx[i]), args.aggregation))	
#             full.apparent <- rbind(full.apparent, full.apparent.i)
            noinf.error.i <- do.call("aggregation.fun", c(list(full.data=actual.data, type="noinf", 
                  response=response , x=x, model=fullmodel, cplx=cplx[i]), args.aggregation))	
            noinf.error <- rbind(noinf.error, noinf.error.i)	
         }
      }
   }

   dimnames(noinf.error) <- NULL
#    dimnames(full.apparent) <- NULL

#    attr(full.apparent, "addattr") <- attributes(full.apparent.i)$addattr
#    fullsample.attr <- attr(full.apparent, "addattr")
#    attr(noinf.error, "addattr") <- attributes(full.apparent.i)$addattr

   pP. <- ls(pattern="predictProb.", name=".GlobalEnv")
   try(eval(call("sfExport", pP., debug=debug)), silent=!debug)
   PLL. <- ls(pattern="PLL.", name=".GlobalEnv")
   try(eval(call("sfExport", PLL., debug=debug)), silent=!debug)
   
   if (environmentName(environment(fit.fun))!="peperr"){sfExport("fit.fun")}
   if (environmentName(environment(aggregation.fun))!="peperr"){sfExport("aggregation.fun")}
   if (environmentName(environment(complexity))!="peperr"){sfExport("complexity")}
   sfExport("response", "x", "sample.index.full", "sample.n", "not.in.sample.full", 
      "args.fit", "args.complexity", "binary")

   if(trace){message("Evaluation on slaves starts now")}

   sample.fun <- function(actual.sample){
      if (trace && actual.sample<sample.n){
         cat("Sample run", actual.sample, "of", (sample.n-1), "\n")
      }
      if (!is.Surv(response)&&!is.matrix(response)){
         response <- matrix(response, ncol=1)
      }
      if (is.function(complexity)){
         sample.complexity <- do.call("complexity",
            c(list(response=response[unique(sample.index.full[[actual.sample]]),],
            x=x[unique(sample.index.full[[actual.sample]]),], full.data=actual.data), args.complexity))
      } else {
         if (is.vector(complexity)){
            sample.complexity <- complexity
         } else sample.complexity <- 0
      }

      if (is.Surv(response)){
         lipec.oob <- c()
         pll.oob <- c()
      }
      actual.error <- c()
      if (is.list(sample.complexity)){
         for (i in 1:length(sample.complexity[[1]])){
            list.sample.complexity <- lapply(sample.complexity, function(arg) arg[i])
            sample.fit <- do.call("fit.fun", c(list(response=response[unique(sample.index.full[[actual.sample]]),],
               x=x[unique(sample.index.full[[actual.sample]]),], 
               cplx=list.sample.complexity), args.fit))
         #class(sample.fit) <- c(class(sample.fit), "peperrinterinternal")

            actual.error.i <-  do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",
               response=response[not.in.sample.full[[actual.sample]],], 
               x=x[not.in.sample.full[[actual.sample]],], model=sample.fit, 
               cplx=list.sample.complexity, fullsample.attr=fullsample.attr), 
               args.aggregation))
            actual.error <- rbind(actual.error, actual.error.i)
 
            if (is.Surv(response)){
               km.fit <- survival::survfit(Surv(time, status)~1,
                  data=actual.data[not.in.sample.full[[actual.sample]],])
               km.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
                  response=response, x=x, model=km.fit), args.aggregation))
               km.pred <- summary(object=km.fit, times=fullsample.attr)$surv 
               km.weight <- -1*diff(km.pred)
               lipec.oob.i <- sum(actual.error.i[1:(length(km.weight))]*km.weight, na.rm=TRUE)
               lipec.oob <- rbind(lipec.oob, lipec.oob.i)
               if (exists(paste("PLL.", class(sample.fit), sep=""))){
                  pll.oob.i <- PLL(object=sample.fit, newdata=x[not.in.sample.full[[actual.sample]],],
                     newtime=time[not.in.sample.full[[actual.sample]]], 
                     newstatus=status[not.in.sample.full[[actual.sample]]], complexity=list.sample.complexity)
                  pll.oob <- rbind(pll.oob, pll.oob.i)
               }
            } else {
               if(binary){
                  logreg.fit <- glm(formula=response~1, data=actual.data[not.in.sample.full[[actual.sample]],],
                     family=binomial())
                  logreg.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",  
                     response=response, x=x, model=logreg.fit), args.aggregation))
                  }
               }
            }
         } else {
            if (is.vector(sample.complexity)){
                for (i in 1:length(sample.complexity)){
                   sample.fit <- do.call("fit.fun",
                      c(list(response=response[unique(sample.index.full[[actual.sample]]),],
                      x=x[unique(sample.index.full[[actual.sample]]),], 
                     cplx=sample.complexity[i]), args.fit))
         #class(sample.fit) <- c(class(sample.fit), "peperrinterinternal")

                   actual.error.i <-  do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",
                      response=response[not.in.sample.full[[actual.sample]],], 
                      x=x[not.in.sample.full[[actual.sample]],], model=sample.fit, 
                      cplx=sample.complexity[i], fullsample.attr=fullsample.attr), 
                      args.aggregation))
                   actual.error <- rbind(actual.error, actual.error.i)

                   if (is.Surv(response)){
                     km.fit <- survival::survfit(Surv(time, status)~1,
                        data=actual.data[not.in.sample.full[[actual.sample]],])
                     km.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent", 
                        response=response, x=x, model=km.fit), args.aggregation))
                     km.pred <- summary(object=km.fit, times=fullsample.attr)$surv 
                     km.weight <- -1*diff(km.pred)
                     lipec.oob.i <- sum(actual.error.i[1:(length(km.weight))]*km.weight, na.rm=TRUE)
                     lipec.oob <- rbind(lipec.oob, lipec.oob.i)
                     if (exists(paste("PLL.", class(sample.fit), sep=""))){
                        pll.oob.i <- PLL(object=sample.fit, newdata=x[not.in.sample.full[[actual.sample]],],
                           newtime=time[not.in.sample.full[[actual.sample]]], 
                           newstatus=status[not.in.sample.full[[actual.sample]]], complexity=sample.complexity[i])
                        pll.oob <- rbind(pll.oob, pll.oob.i)
                     }
                  } else {
                     if (binary){
                     logreg.fit <- glm(formula=response~1, data=actual.data[not.in.sample.full[[actual.sample]],],
                        family=binomial())
                     logreg.apparent <- do.call("aggregation.fun", c(list(full.data=actual.data, type="apparent",  
                        response=response, x=x, model=logreg.fit), args.aggregation))
                     }
                     }
                  }
              }
          }
      dimnames(actual.error) <- NULL
      
      if (is.Surv(response)){
         dimnames(lipec.oob) <- NULL
         dimnames(pll.oob) <- NULL
         out <- list(actual.error=actual.error, sample.complexity=sample.complexity, 
            sample.fit=sample.fit,
            lipec.oob=lipec.oob, pll.oob=pll.oob, km.apparent=km.apparent)
      } else {
         if(binary){
            out <- list(actual.error=actual.error, sample.complexity=sample.complexity,
                      sample.fit=sample.fit, logreg.apparent=logreg.apparent)
         } else {
            out <- list(actual.error=actual.error, sample.complexity=sample.complexity,
               sample.fit=sample.fit)
         }
      }
   out
}

   sample.error.list <- sfClusterApplyLB(as.list(1:sample.n), sample.fun)
   sample.error.full <- lapply(sample.error.list, function(arg) arg$actual.error)
   sample.complexity.full <- unlist(lapply(sample.error.list, function(arg) arg$sample.complexity))
   sample.fit.full <- lapply(sample.error.list, function(arg) arg$sample.fit)
   if (is.Surv(response)){
   sample.lipec.full <- lapply(sample.error.list, function(arg) arg$lipec.oob)
   sample.pll.full <- lapply(sample.error.list, function(arg) arg$pll.oob)
   sample.km.full <- lapply(sample.error.list, function(arg) arg$km.apparent)
   attr(null.model, "addattr") <- NULL
   } else {
   if (binary){
   sample.lrm.full <- lapply(sample.error.list, function(arg) arg$logreg.fit)
   sample.null.model.full <- lapply(sample.error.list, function(arg) arg$logreg.apparent)
   }
   }

   sfStop(nostop=noclusterstop)

   if (is.vector(sample.complexity.full)){
      sample.complexity <- sample.complexity.full[-sample.n]
   }
   if (is.list(cplx)){
      sample.complexity <- sample.complexity.full[1:(length(sample.complexity.full)-length(cplx))]
   }

   if (!is.function(complexity) && !is.list(complexity)){
   sample.complexity <- unique(sample.complexity)
   }


   if (is.Surv(response)){
   output <- list(indices=list(sample.index=sample.index, not.in.sample=not.in.sample),
      selected.complexity=cplx, complexity=complexity, 
      response=response, full.model.fit=sample.fit.full[[sample.n]],
      full.apparent=sample.error.full[[sample.n]], noinf.error=noinf.error, 
      attribute=fullsample.attr,
      sample.error=sample.error.full[1:(sample.n-1)], sample.complexity=sample.complexity,  
      sample.lipec=sample.lipec.full[1:(sample.n-1)], sample.pll=sample.pll.full[1:(sample.n-1)],
      null.model=null.model, null.model.fit=km.fit, sample.null.model=sample.km.full[1:(sample.n-1)])
   } else {
      if (binary){
         output <- list(indices=list(sample.index=sample.index, not.in.sample=not.in.sample),
            selected.complexity=cplx, complexity=complexity, 
            response=response, full.model.fit=sample.fit.full[[sample.n]],
            full.apparent=sample.error.full[[sample.n]], noinf.error=noinf.error,
            attribute=fullsample.attr,
            sample.error=sample.error.full[1:(sample.n-1)], sample.complexity=sample.complexity,
            null.model=null.model, null.model.fit=logreg.fit, sample.null.model=sample.null.model.full[1:(sample.n-1)])
      } else {
          output <- list(indices=list(sample.index=sample.index, not.in.sample=not.in.sample),
            selected.complexity=cplx, complexity=complexity, 
            response=response, full.model.fit=sample.fit.full[[sample.n]],
            full.apparent=sample.error.full[[sample.n]], noinf.error=noinf.error,
            attribute=fullsample.attr,
            sample.error=sample.error.full[1:(sample.n-1)], sample.complexity=sample.complexity)
      }
   }
   class(output) <- "peperr"
   output
}


