# partial_bdd gives the estimates of conditinal probablity
partial_bdd <- function(cond_prob_array, bdd_type=c("B-P", "Siddique")){
  p <- q <- array(0, c(length(cond_prob_array[1,1,1,]), 8))
  if (bdd_type=="B-P"){
    p[,1] <- cond_prob_array[1,1,1,]+cond_prob_array[2,2,2,]-1
    p[,2] <- cond_prob_array[1,1,2,]+cond_prob_array[2,2,2,]-1
    p[,3] <- cond_prob_array[2,2,1,]+cond_prob_array[1,1,2,]-1
    p[,4] <- cond_prob_array[1,1,1,]+cond_prob_array[2,2,1,]-1
    p[,5] <- 2*cond_prob_array[1,1,1,]+cond_prob_array[2,2,1,]+cond_prob_array[2,1,2,]+cond_prob_array[2,2,2,]-2
    p[,6] <- cond_prob_array[1,1,1,]+2*cond_prob_array[2,2,1,]+cond_prob_array[1,1,2,]+cond_prob_array[1,2,2,]-2
    p[,7] <- cond_prob_array[2,1,1,]+cond_prob_array[2,2,1,]+2*cond_prob_array[1,1,2,]+cond_prob_array[2,2,2,]-2
    p[,8] <- cond_prob_array[1,1,1,]+cond_prob_array[1,2,1,]+cond_prob_array[1,1,2,]+2*cond_prob_array[2,2,2,]-2
    
    q[,1] <- 1-cond_prob_array[2,1,1,]-cond_prob_array[1,2,2,]
    q[,2] <- 1-cond_prob_array[1,2,1,]-cond_prob_array[2,1,2,]
    q[,3] <- 1-cond_prob_array[1,2,1,]-cond_prob_array[2,1,1,]
    q[,4] <- 1-cond_prob_array[1,2,2,]-cond_prob_array[2,1,2,]
    q[,5] <- 2-2*cond_prob_array[1,2,1,]-cond_prob_array[2,1,1,]-cond_prob_array[2,1,2,]-cond_prob_array[2,2,2,]
    q[,6] <- 2-cond_prob_array[1,2,1,]-2*cond_prob_array[2,1,1,]-cond_prob_array[1,1,2,]-cond_prob_array[1,2,2,]
    q[,7] <- 2-cond_prob_array[2,1,1,]-cond_prob_array[2,2,1,]-2*cond_prob_array[1,2,2,]-cond_prob_array[2,1,2,]
    q[,8] <- 2-cond_prob_array[1,1,1,]-cond_prob_array[1,2,1,]-cond_prob_array[1,2,2,]-2*cond_prob_array[2,1,2,]
  }
  list(p=p, q=q)
}

# get_cond_prob returns the estimated conditional probability matrix given the training samples
get_cond_prob <- function(y, trt, z, x, newx=x){
  cond_prob <- array(0, c(2,2,2,NROW(newx)))
  for (y_index in c(1,2)){
    for (trt_index in c(1,2)){
      for (z_index in c(1,2)){
        y_transformed <- (y==2*(y_index-3/2)) * (trt==2*(trt_index-3/2))
        cond_prob[y_index, trt_index, z_index,] <- ks(xx=x[z==2*(z_index-3/2),], yy=y_transformed[z==2*(z_index-3/2)], xx.test=newx)
      }
    }
  }
  cond_prob
}

# get_cond_prob returns the estimated conditional probability matrix given the training samples
get_z_prob <- function(z, x, newx=x){
  cond_prob <- ks(xx=x, yy=(z+1)/2, xx.test=newx)
  cond_prob
}

# get_debiased_vector returns the vectors to construct the debiased term
get_vector_to_debias <- function(term=c("p", "q"), number=c(1:8)){
  if (term=="p"){
    if (number==1){
      c(1,0,0,0,0,0,0,1)
    } else if (number==2){
      c(0,1,0,0,0,0,0,1)
    } else if (number==3){
      c(0,1,0,0,0,0,1,0)
    } else if (number==4){
      c(1,0,0,0,0,0,1,0)
    } else if (number==5){
      c(2,0,0,0,0,1,1,1)
    } else if (number==6){
      c(1,1,0,1,0,0,2,0)
    } else if (number==7){
      c(0,2,0,0,1,0,1,1)
    } else {
      c(1,1,1,0,0,0,0,2)
    }
  } else {
    if (number==1){
      c(0,0,0,-1,-1,0,0,0)
    } else if (number==2){
      c(0,0,-1,0,0,-1,0,0)
    } else if (number==3){
      c(0,0,-1,0,-1,0,0,0)
    } else if (number==4){
      c(0,0,0,-1,0,-1,0,0)
    } else if (number==5){
      c(0,0,-2,0,-1,-1,0,-1)
    } else if (number==6){
      c(0,-1,-1,-1,-2,0,0,0)
    } else if (number==7){
      c(0,0,0,-2,-1,-1,-1,0)
    } else {
      c(-1,0,-1,-1,0,-2,0,0)
    }
  }
}

# get_debiased_vector returns the vectors to construct the debiased term
get_y_to_debias <- function(y_a=c(1,1)){
 if (y_a[1] == -1 & y_a[2] == -1){
   c(1,1,0,0,0,0,0,0)
 } else if (y_a[1] == -1 & y_a[2] == 1){
   c(0,0,1,1,0,0,0,0)
 } else if (y_a[1] == 1 & y_a[2] == -1){
   c(0,0,0,0,1,1,0,0)
 } else {
   c(0,0,0,0,0,0,1,1)
 }
}

# get_debiased_vector returns the vectors to construct the debiased term
get_z_to_debias <- function(z=1){
  if (z==-1){
    c(1,0,1,0,1,0,1,0)
  } else {
    c(0,1,0,1,0,1,0,1)
  }
}
