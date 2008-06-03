`extract.fun` <-
function(funs=NULL) 
{
   if (is.null(funs)){
      extract.packages <- (.packages())
      functions <- c()
      for (i in seq(along=ls(name=".GlobalEnv"))){
         name <- try(get(as.character(ls(name=".GlobalEnv")[i])), silent=TRUE)
         packagename <- environmentName(environment(name))
         if (packagename=="R_GlobalEnv"){
            functions <- c(functions, ls(name=".GlobalEnv")[i])
         }
      }
   } else {
      packages <- c()
      functions <- c()
      for (j in seq(along=funs)){
         if (is.function(funs[[j]])){
            packages <- c(packages, extract(funs[[j]])$packages)
            functions <- c(functions, extract(funs[[j]])$functions)
         }
      }
   }
   extract.packages <- unique(packages)
   if (length(extract.packages)==0) extract.packages <- NULL

   extract.functions <- unique(functions)
   
   extract.fun.list <- list(packages=extract.packages, functions=extract.functions)
   extract.fun.list
}

