cf.iv_itr <- function(y, trt, z, x, split_index = NULL, lambda_seq=10^(seq(-2,0, length.out=5)), 
                      bandwidth_seq=10^(seq(-2,1, length.out=10)), d=0.5 , bdd_type=c("B-P", "Siddique"), 
                      max_iter=10, debias=TRUE, parallel=FALSE, solver="CVXR", nfold=5){
  sample_size <- length(y)
  if (is.null(split_index)){
    split_index <- sample(1:nfold, sample_size, replace=TRUE)
  }
  
  if (!parallel){
    cf.error <- array(0, c(length(lambda_seq), length(bandwidth_seq)))
    for (j in 1:length(bandwidth_seq)){
      for (i in 1:length(lambda_seq)){
        loss <- sapply(1:nfold, function(k){
          tryCatch({
            fit <- iv_itr(y=y, trt=trt, z=z, x=x, split_index = (split_index!=k), lambda=lambda_seq[i], d=d,
                          bandwidth=bandwidth_seq[j], bdd_type=bdd_type, max_iter=max_iter, debias=debias,
                          solver = solver)
            value_bdd(fit=fit,y=y, trt=trt, z=z, x=x, split_index = (split_index!=k), d=d)
          }, error=function(e){
            warning("Skipping failed CV fold: lambda=", lambda_seq[i],
                    ", bandwidth=", bandwidth_seq[j], ", fold=", k,
                    " | ", conditionMessage(e), call. = FALSE)
            NA_real_
          })
        })
        cf.error[i,j] <- if (all(is.na(loss))) Inf else mean(loss, na.rm = TRUE)
      }
    }
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE) #as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE',unset='1'))
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    cf.error <- foreach(j=1:length(bandwidth_seq), .errorhandling = "remove",.combine = rbind)%dopar%{
      cf.error <- array(0, c(1, length(lambda_seq)))
      for (i in 1:length(lambda_seq)){
        loss <- sapply(1:nfold, function(k){
          fit <- iv_itr(y=y, trt=trt, z=z, x=x, split_index = (split_index!=k), lambda=lambda_seq[i], d=d, 
                        bandwidth=bandwidth_seq[j], bdd_type=bdd_type, max_iter=max_iter, debias=debias,
                        solver = solver)
          value_bdd(fit=fit,y=y, trt=trt, z=z, x=x, split_index = (split_index!=k), d=d)
        })
        cf.error[1,i] <- mean(loss, na.rm = TRUE)
      }
      cf.error
    }
    stopCluster(cl)
    cf.error <- t(cf.error)
  }
  
  finite.index <- which(is.finite(cf.error), arr.ind = TRUE)
  if (nrow(finite.index) == 0){
    stop("All cross-validation fits failed.", call. = FALSE)
  }
  finite.error <- cf.error[finite.index]
  finite.index <- finite.index[order(finite.error),, drop=FALSE]
  min.index <- finite.index[1,, drop=FALSE]
  if (!parallel){
    fit <- NULL
    last_error <- NULL
    for (candidate_index in seq_len(nrow(finite.index))){
      min.index <- finite.index[candidate_index,, drop=FALSE]
      fit <- tryCatch({
        lapply(1:nfold, function(k){
          fit <- iv_itr(y=y, trt=trt, z=z, x=x, split_index = (split_index!=k),
                 lambda=lambda_seq[min.index[1]], d=d, bandwidth=bandwidth_seq[min.index[2]],
                 bdd_type=bdd_type, max_iter=max_iter, debias=debias,solver = solver)
          fit
        })
      }, error=function(e){
        last_error <<- e
        NULL
      })
      if (!is.null(fit)){
        break
      }
    }
    if (is.null(fit)){
      stop("All final fits failed after cross-validation. Last error: ",
           conditionMessage(last_error), call. = FALSE)
    }
  } else {
      library(doParallel)
      n_cores <- detectCores(all.tests = FALSE, logical = TRUE) #as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE',unset='1'))
      cl <- makeCluster(n_cores)
      registerDoParallel(cl)
      fit<- foreach(k=1:nfold, .errorhandling = "remove")%dopar%{
        fit <- iv_itr(y=y, trt=trt, z=z, x=x, split_index = (split_index!=k), lambda=lambda_seq[min.index[1]], d=d, 
               bandwidth=bandwidth_seq[min.index[2]], bdd_type=bdd_type, max_iter=max_iter, debias=debias,
               solver = solver)
        fit
      }
      stopCluster(cl)
  }
  
  list(fit=fit, cv.error=cf.error, bandwidth_opt=bandwidth_seq[min.index[2]], lambda_opt=lambda_seq[min.index[1]], bdd_type=bdd_type, nfold=nfold)
}
