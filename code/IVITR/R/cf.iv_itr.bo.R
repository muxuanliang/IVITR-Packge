cf.iv_itr.bo <- function(y, trt, z, x, split_index = NULL, lambda_seq=10^(seq(-2,0, length.out=5)),
                      bandwidth_seq=10^(seq(-2,1, length.out=10)), bdd_type=c("B-P", "Siddique")){
  sample_size <- length(y)
  if (is.null(split_index)){
    split_index <- rnorm(sample_size)>0
  }

  cf.error <- array(0, c(length(lambda_seq), length(bandwidth_seq)))
  for (i in 1:length(lambda_seq)){
    for (j in 1:length(bandwidth_seq)){
      fit1 <- iv_itr.bo(y=y, trt=trt, z=z, x=x, split_index = split_index, lambda=lambda_seq[i],
                     bandwidth=bandwidth_seq[j], bdd_type=bdd_type)
      loss1 <- value_bdd(fit=fit1,y=y, trt=trt, z=z, x=x, split_index = split_index, d=1/2)
      fit2 <- iv_itr.bo(y=y, trt=trt, z=z, x=x, split_index = !split_index, lambda=lambda_seq[i],
                     bandwidth=bandwidth_seq[j], bdd_type=bdd_type)
      loss2 <- value_bdd(fit=fit2,y=y, trt=trt, z=z, x=x, split_index = !split_index, d=1/2)
      cf.error[i,j] <- (loss1+loss2)/2
    }
  }
  min.index <- which(cf.error==min(cf.error), arr.ind = TRUE)
  fit <- NULL
  fit[[1]] <- iv_itr.bo(y=y, trt=trt, z=z, x=x, split_index = split_index, lambda=lambda_seq[min.index[1]],
                     bandwidth=bandwidth_seq[min.index[2]], bdd_type=bdd_type)
  fit[[2]] <- iv_itr.bo(y=y, trt=trt, z=z, x=x, split_index = !split_index, lambda=lambda_seq[min.index[1]],
                     bandwidth=bandwidth_seq[min.index[2]], bdd_type=bdd_type)
  list(fit=fit, cv.error=cf.error, bandwidth_opt=bandwidth_seq[min.index[2]], lambda_opt=lambda_seq[min.index[1]], bdd_type=bdd_type)
}
