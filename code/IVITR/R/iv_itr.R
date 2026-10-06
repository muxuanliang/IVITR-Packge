iv_itr <- function(y, trt, z, x, split_index = NULL, lambda=0.01, d=0.5, bandwidth=0.1, 
                   bdd_type=c("B-P", "Siddique"), max_iter=10, max_iter_inner=1000, debias=TRUE, solver="CVXR"){
  sample_size <- length(y)
  if (is.null(split_index)){
    split_index <- rnorm(sample_size)>0
  }
  data_ns <- list(y=y[split_index], trt=trt[split_index], z=z[split_index], x=x[split_index,])
  data_solve <- list(y=y[!split_index], trt=trt[!split_index], z=z[!split_index], x=x[!split_index,])

  cond_prob <- get_cond_prob(y=data_ns$y, trt=data_ns$trt, z=data_ns$z, x=data_ns$x, newx=data_solve$x)
  pq <- partial_bdd(cond_prob_array = cond_prob, bdd_type = bdd_type)
  k_matrix <- rbf_gauss(data_solve$x , bandwidth = bandwidth)
  p_max <- apply(pq$p,1,max)
  q_min <- apply(pq$q,1,min)
  positive_weights <- abs(q_min)*(q_min>=0)
  negative_weights <- abs(p_max)*(p_max<=0)

  # which max or min
  debiased_z_weight <- t(sapply(data_solve$z, function(t){get_z_to_debias(t)}))
  cond_prob_z <- get_z_prob(z=data_ns$z, x=data_ns$x, newx=data_solve$x)
  debiased_y <- (cbind(cond_prob[1,1,1,], cond_prob[1,1,2,], cond_prob[1,2,1,],cond_prob[1,2,2,],
                       cond_prob[2,1,1,], cond_prob[2,1,2,], cond_prob[2,2,1,],cond_prob[2,2,2,])-
                   t(apply(cbind(data_solve$y,data_solve$trt), 1, function(t){get_y_to_debias(t)}))) * debiased_z_weight

  which_q_min <- apply(pq$q,1,which.min)
  vector_to_debias_q <- t(sapply(which_q_min, function(t){get_vector_to_debias('q',t)}))
  positive_debiased <- (q_min>=0) * apply(vector_to_debias_q * debiased_y,1,sum) * (0.5+(cond_prob_z-0.5) * data_solve$z)^(-1)
  which_p_max <- apply(pq$p,1,which.max)
  vector_to_debias_p <- t(sapply(which_p_max, function(t){get_vector_to_debias('p',t)}))
  negative_debiased <- (-1)*(p_max<=0) * apply(vector_to_debias_p * debiased_y,1,sum) * (0.5+(cond_prob_z-0.5) * data_solve$z)^(-1)

  positive_weights <- abs(q_min)*(q_min>=0)-positive_debiased
  negative_weights <- abs(p_max)*(p_max<=0)-negative_debiased

  # warm-up
  fit_init <- iv_itr_warm(y, trt, z, x, split_index = split_index, lambda=lambda, 
                          d=d, bandwidth=bandwidth, bdd_type=bdd_type, solver = solver)
  alpha_updated <- fit_init$alpha_est
  iter <- 0
  diff_value <- 1
  
  if(!debias){
    return(list(alpha_est=alpha_updated * (abs(alpha_updated)>1e-5), basis_x=data_solve$x, 
         status="optimal", alpha_warm=alpha_updated * (abs(alpha_updated)>1e-5), bandwidth=bandwidth,
         bdd_type=bdd_type))
  }
  
  if (solver=="HD"){
    alpha_positive_prev1 <- fit_init$alpha_positive_prev1
    alpha_positive_prev2 <- fit_init$alpha_positive_prev2
    alpha_negative_prev1 <- fit_init$alpha_negative_prev1
    alpha_negative_prev2 <- fit_init$alpha_negative_prev2
  }
  
  # use DCA algorithm
  positive_weights_g <- positive_weights*(positive_weights>0)
  negative_weights_g <- negative_weights*(negative_weights>0)
  positive_weights_h <- -positive_weights*(positive_weights<0)
  negative_weights_h <- -negative_weights*(negative_weights<0)
  while(iter <= max_iter & abs(diff_value)>1e-5){
    alpha_prev <- alpha_updated
    ## P-problem
   
    
    if (solver=="CVXR"){
      yt <- k_matrix %*% (positive_weights_h*(k_matrix %*%alpha_prev<=1)*(k_matrix %*%alpha_prev>0)*(-1)+
                            positive_weights_h*(k_matrix %*%alpha_prev<=0)*(1-1/d)-
                            negative_weights_h*(k_matrix %*%alpha_prev>= -1)*(k_matrix %*%alpha_prev<0)*(-1)-
                            negative_weights_h*(k_matrix %*%alpha_prev>=0)*(1-1/d))/NROW(data_solve$x)
      # square-root of K
      svd_k <- svd(k_matrix)
      chol_K <- diag(sqrt(svd_k$d)+0.01) %*% t(svd_k$u)
      ## Q-problem
      alpha <- CVXR::Variable(sum(!split_index))
      slack_positive <- CVXR::Variable(sum(!split_index))
      slack_negative <- CVXR::Variable(sum(!split_index))
      objective <- CVXR::Minimize(mean(positive_weights_g*slack_positive)+
                                    mean(negative_weights_g*slack_negative)+
                                    lambda* sum((chol_K %*% alpha)^2)-sum(alpha* yt))
      constraint1 <- slack_positive >= 0
      constraint2 <- slack_positive >= (1-k_matrix %*% alpha)
      constraint3 <- slack_positive >= (1-(1-d)/d * k_matrix %*% alpha)
      constraint4 <- slack_negative >= 0
      constraint5 <- slack_negative >= (1+k_matrix %*% alpha)
      constraint6 <- slack_negative >= (1+(1-d)/d * k_matrix %*% alpha)
      problem <- CVXR::Problem(objective, constraints = list(constraint1, constraint2,constraint3, constraint4,
                                                             constraint5, constraint6))
      result <- CVXR::solve(problem, solver = "OSQP")
      alpha_updated <- cvxr_value(result, alpha)
      
      diff_value <- mean(positive_weights*pmax(0,(1-k_matrix %*% alpha_updated ),(1-(1-d)/d * k_matrix %*% alpha_updated)))+
        mean(negative_weights*pmax(0,(1+k_matrix %*% alpha_updated),(1+(1-d)/d * k_matrix %*% alpha_updated)))+
        lambda* t(alpha_updated) %*% k_matrix %*% alpha_updated-mean(positive_weights*pmax(0,(1-k_matrix %*% alpha_prev),(1-(1-d)/d * k_matrix %*% alpha_prev)))-
        mean(negative_weights*pmax(0,(1+k_matrix %*% alpha_prev),(1+(1-d)/d * k_matrix %*% alpha_prev)))-
        lambda* t(alpha_prev) %*% k_matrix %*% alpha_prev
      
      iter <- iter+1
    } else if (solver=="HD"){
      # square-root of K
      yt <- (positive_weights_h*(k_matrix %*%alpha_prev<=1)*(k_matrix %*%alpha_prev>0)*(-1)+
                            positive_weights_h*(k_matrix %*%alpha_prev<=0)*(1-1/d)-
                            negative_weights_h*(k_matrix %*%alpha_prev>= -1)*(k_matrix %*%alpha_prev<0)*(-1)-
                            negative_weights_h*(k_matrix %*%alpha_prev>=0)*(1-1/d))/NROW(data_solve$x)
      iter_inner <- 0
      diff_value_inner <- 1
      batch <- floor(sum(!split_index)/1000)+1
      sample_index <- sample(1:batch, sum(!split_index), replace = TRUE)
      alpha_prev_inner <- alpha_prev
      while(iter_inner <= max_iter_inner & abs(diff_value_inner)>1e-5){
        update_index <- iter_inner%%batch+1
        svd_k <- svd(k_matrix[sample_index==update_index, sample_index==update_index])
        chol_K <- diag(sqrt(svd_k$d)+0.01) %*% t(svd_k$u)
        # use convex optimization
        alpha_positive1 <- CVXR::Variable(sum(sample_index==update_index))
        alpha_negative1 <- CVXR::Variable(sum(sample_index==update_index))
        alpha_positive2 <- CVXR::Variable(sum(sample_index==update_index))
        alpha_negative2 <- CVXR::Variable(sum(sample_index==update_index))
        alpha <- CVXR::Variable(sum(sample_index==update_index))
        objective <- CVXR::Minimize(0.5* sum((chol_K %*% alpha)^2)
                                    +t(alpha_positive_prev1[sample_index!=update_index]+(1-d)/d *alpha_positive_prev2[sample_index!=update_index] -alpha_negative_prev1[sample_index!=update_index]-(1-d)/d *alpha_negative_prev2[sample_index!=update_index]+yt[sample_index!=update_index]/(2*lambda)) %*% 
                                      k_matrix[sample_index!=update_index, sample_index==update_index]%*% alpha-
                                      sum(alpha_positive1+alpha_negative1+alpha_positive2+alpha_negative2))
        constraint0 <- alpha==alpha_positive1+(1-d)/d*alpha_positive2-alpha_negative1-(1-d)/d *alpha_negative2+yt[sample_index==update_index]/(2*lambda)
        constraint1 <- alpha_positive1 >= 0
        constraint2 <- alpha_positive2 >= 0
        constraint3 <- (alpha_positive1+alpha_positive2) <= positive_weights_g[sample_index==update_index]/(2*lambda*sum(!split_index))
        constraint4 <- alpha_negative1 >= 0
        constraint5 <- alpha_negative2 >= 0
        constraint6 <- (alpha_negative1+alpha_negative2) <= negative_weights_g[sample_index==update_index]/(2*lambda*sum(!split_index))
        problem <- CVXR::Problem(objective, constraints = list(constraint0,constraint1, constraint2,constraint3, constraint4,constraint5, constraint6))
        result <- CVXR::solve(problem, solver = "OSQP")
        alpha_positive_updated1 <- alpha_positive_prev1
        alpha_negative_updated1 <- alpha_negative_prev1
        alpha_positive_updated2 <- alpha_positive_prev2
        alpha_negative_updated2 <- alpha_negative_prev2
        alpha_positive_updated1[sample_index==update_index] <- cvxr_value(result, alpha_positive1)
        alpha_negative_updated1[sample_index==update_index] <- cvxr_value(result, alpha_negative1)
        alpha_positive_updated2[sample_index==update_index] <- cvxr_value(result, alpha_positive2)
        alpha_negative_updated2[sample_index==update_index] <- cvxr_value(result, alpha_negative2)
        alpha_updated_inner <- alpha_prev_inner
        alpha_updated_inner[sample_index==update_index] <- cvxr_value(result, alpha)
        diff_value_inner <- mean(positive_weights_g*pmax(0,(1-k_matrix %*% alpha_updated_inner),(1-(1-d)/d * k_matrix %*% alpha_updated_inner)))+
          mean(negative_weights_g*pmax(0,(1+k_matrix %*% alpha_updated_inner),(1+(1-d)/d * k_matrix %*% alpha_updated_inner)))+
          lambda* t(alpha_updated_inner) %*% k_matrix %*% alpha_updated_inner-sum(alpha_updated_inner* k_matrix %*% yt)-
        mean(positive_weights_g*pmax(0,(1-k_matrix %*% alpha_prev_inner),(1-(1-d)/d * k_matrix %*% alpha_prev_inner)))-
          mean(negative_weights_g*pmax(0,(1+k_matrix %*% alpha_prev_inner),(1+(1-d)/d * k_matrix %*% alpha_prev_inner)))-
          lambda* t(alpha_prev_inner) %*% k_matrix %*% alpha_prev_inner+sum(alpha_prev_inner* k_matrix %*% yt)
        if (diff_value_inner<0){
          alpha_prev_inner <- alpha_updated_inner
          alpha_positive_prev1 <- alpha_positive_updated1
          alpha_positive_prev2 <- alpha_positive_updated2
          alpha_negative_prev1 <- alpha_negative_updated1
          alpha_negative_prev2 <- alpha_negative_updated2
          iter_inner <- iter_inner+1
        } else {
          break
        }
      }
      alpha_updated <- alpha_prev_inner
      diff_value <- mean(positive_weights*pmax(0,(1-k_matrix %*% alpha_updated ),(1-(1-d)/d * k_matrix %*% alpha_updated)))+
        mean(negative_weights*pmax(0,(1+k_matrix %*% alpha_updated),(1+(1-d)/d * k_matrix %*% alpha_updated)))+
        lambda* t(alpha_updated) %*% k_matrix %*% alpha_updated-mean(positive_weights*pmax(0,(1-k_matrix %*% alpha_prev),(1-(1-d)/d * k_matrix %*% alpha_prev)))-
        mean(negative_weights*pmax(0,(1+k_matrix %*% alpha_prev),(1+(1-d)/d * k_matrix %*% alpha_prev)))-
        lambda* t(alpha_prev) %*% k_matrix %*% alpha_prev
      
      iter <- iter+1
    }
  }
  if (iter > max_iter){
    status <- "not converge"
  } else {
    status <- "optimal"
  }
  
  list(alpha_est=alpha_updated * (abs(alpha_updated)>1e-5), basis_x=data_solve$x, 
       status=status, alpha_warm=fit_init$alpha_est, bandwidth=bandwidth,
       bdd_type=bdd_type)
}
