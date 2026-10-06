predict.iv_itr <- function(fit, newx, type=c("link", "decision"), d=0.5){
  predict.iv_itr_split <- sapply(1:fit$nfold, function(k){
    k_matrix <- rbf_gauss(x=fit$fit[[k]]$basis_x, newx=newx, bandwidth = fit$bandwidth_opt)
    k_matrix %*% fit$fit[[k]]$alpha_est
  })
  predict.iv_itr <- apply(predict.iv_itr_split, 1, mean)
  if (type=="link"){
    return(predict.iv_itr)
  } else {
    if (d==0.5){
      return((predict.iv_itr>0)-(predict.iv_itr< 0))
    } else {
      return((predict.iv_itr> (0.5-d))-(predict.iv_itr< -(0.5-d)))
    }
  }
}