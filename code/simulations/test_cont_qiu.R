sim <- function(n=500, eta=0.5, setting1=1, setting2=1, d=0.5){
  library(doParallel)
  n_cores <- detectCores(all.tests = FALSE, logical = TRUE) #as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE',unset='1'))
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  p <- 5
  lambda <- xi <- delta <- eta
  res_proposed <- read_csv(paste0(getwd(),"/simulation/res_",n,"_",lambda,"_", xi, "_", delta, "_",setting1, "_", setting2, "_", d,"_cont.csv"))
  
  print(system.time(res <- foreach(i=1:200, .packages="IVITR", .errorhandling = "remove",.combine = rbind)%dopar%{
    set.seed(i)
    
    x <- array(runif(n*p, min=-1, max=1), c(n, p))
    u <- runif(n, min=-1, max=1)
    z <- (rnorm(n)>0)
    x_u_z_cmd <- cbind(x, u, z)
    trt <- apply(x_u_z_cmd,1, function(t){rbinom(1,1,prob = exp(8*t[p+2]+t[1]-7*t[2]+lambda*(1+t[1])*t[p+1])/(1+exp(8*t[p+2]+t[1]-7*t[2]+lambda*(1+t[1])*t[p+1])))})
    x_u_trt_cmd <- cbind(x, u, trt)
    
    get_y_sim <- function(t){
      if (setting1 ==1){
        g1 <- 1-t[1]-t[2]+xi*t[p+1]
      } else {
        g1 <- 1-(t[1])^2+(t[2])^2+xi*t[1]*t[2]*t[p+1]
      }
      if (setting2 ==1){
        g2 <- 0.442*(1-t[1]-t[2]+delta*t[p+1])*t[p+2]
      } else {
        g2 <- (t[2]-delta*0.25*(t[1])^2+4*(t[2]>0)-6*(t[2]<=0))*t[p+2]
      }
      noise <- 0.1*rnorm(1,0,1)
      exp(0.5*g1+0.5*g2+noise+1)/(1+exp(0.5*g1+0.5*g2+noise+1))
    }
    y <- apply(x_u_trt_cmd,1, get_y_sim)
    
    # flip 0,1 to -1, 1
    z_tfd <- 2*(z-0.5)
    trt_tfd <- 2*(trt-0.5)
    y_tfd <- 2*(y-0.5)
    
    # fit
    #fit <- cf.iv_itr(y_tfd, trt_tfd,z_tfd,x, split_index = NULL, d=d, bdd_type="B-P")
    
    # new data
    newx <- array(runif(1000*p, min=-1, max=1), c(1000, p))
    newu <- runif(1000, min=-1, max=1)
    get_true_y_sim <- function(t){
      if (setting1 ==1){
        g1 <- 1-t[1]-t[2]+xi*t[p+1]
      } else {
        g1 <- 1-(t[1])^2+(t[2])^2+xi*t[1]*t[2]*t[p+1]
      }
      if (setting2 ==1){
        g2 <- 0.442*(1-t[1]-t[2]+delta*t[p+1])*t[p+2]
      } else {
        g2 <- (t[2]-delta*0.25*(t[1])^2+4*(t[2]>0)-6*(t[2]<=0))*t[p+2]
      }
      noise <- 0.1*rnorm(1,0,1)
      exp(0.5*g1+0.5*g2+noise+1)/(1+exp(0.5*g1+0.5*g2+noise+1))
    }
    y1 <- apply(x_u_trt_cmd <- cbind(newx, newu, 1),1, get_true_y_sim)
    y0 <- apply(x_u_trt_cmd <- cbind(newx, newu, 0),1, get_true_y_sim)
    
    pr_fit <- res_proposed$value[res_proposed$index==i & res_proposed$metric=="porp rejected" & res_proposed$method=="Debiased"]
    #fit_pred <- predict.iv_itr(fit, newx=newx, type="link")
    fit_qiu <- iv_itr.qiu(y_tfd, trt_tfd,z_tfd,x, newx=newx)
    
    # find cutoff
    cutoff.qiu <- quantile(abs(fit_qiu$delta_eff), probs = pr_fit)
    
    # value on non-rejected
    if (d==0.5){
      value_fit_nonreject <- res_proposed$value[res_proposed$index==i & res_proposed$metric=="value on non-rejected" & res_proposed$method=="Debiased"]
      value_fit_qiu_nonreject <- sum(y1*(fit_qiu$delta_eff>cutoff.qiu)+y0*(fit_qiu$delta_eff< -cutoff.qiu))/sum((fit_qiu$delta_eff>cutoff.qiu)+(fit_qiu$delta_eff<cutoff.qiu))
    } else {
      value_fit_nonreject <- res_proposed$value[res_proposed$index==i & res_proposed$metric=="value on non-rejected" & res_proposed$method=="Debiased"]
      value_fit_qiu_nonreject <- sum(y1*(fit_qiu$delta_eff>cutoff.qiu)+y0*(fit_qiu$delta_eff< -cutoff.qiu))/sum((fit_qiu$delta_eff>cutoff.qiu)+(fit_qiu$delta_eff<cutoff.qiu))
    }
    
    # loss
    newz <- (rnorm(1000)>0)
    x_u_z_cmd_new <- cbind(newx, newu, newz)
    newtrt <- apply(x_u_z_cmd_new,1, function(t){rbinom(1,1,prob = exp(8*t[p+2]+t[1]-7*t[2]+lambda*(1+t[1])*t[p+1])/(1+exp(8*t[p+2]+t[1]-7*t[2]+lambda*(1+t[1])*t[p+1])))})
    x_u_trt_cmd_new <- cbind(newx, newu, newtrt)
    newy <- apply(x_u_trt_cmd_new,1, get_y_sim)
    newv <- runif(1000)
    newy_v <- (newy>= newv)
    newz_tfd <- 2*(newz-0.5)
    newtrt_tfd <- 2*(newtrt-0.5)
    newy_tfd <- 2*(newy_v-0.5)
    
    cond_prob <- get_cond_prob(y=newy_tfd, trt=newtrt_tfd, z=newz_tfd, x=newx)
    prob_bdd <- partial_bdd(cond_prob, bdd_type = "B-P")
    mean(prob_bdd$p*prob_bdd$q>=0)
    mean(prob_bdd$p*prob_bdd$q>=0 & prob_bdd$p<0)
    
    loss_fit <- res_proposed$value[res_proposed$index==i & res_proposed$metric=="loss in value" & res_proposed$method=="Debiased"]
    loss_fit.qiu <- value_bdd.qiu(fit_qiu$delta_eff, cutoff.qiu,newy_tfd, newtrt_tfd, newz_tfd, newx, d)
    
    res <- data.frame(value=c(value_fit_nonreject, value_fit_qiu_nonreject, pr_fit,pr_fit,loss_fit,loss_fit.qiu),
                      d=d, metric=rep(c("value on non-rejected", "porp rejected", "loss in value"), each=2),
                      method=c("Debiased","Qiu"), index=i, n=n, lambda=lambda, xi=xi, delta=delta,
                      setting1=setting1, setting2=setting2)
    res
  }))
  stopCluster(cl)
  write.csv(res, file=paste0(getwd(),"/simulation/res_",n,"_",lambda,"_", xi, "_", delta, "_",setting1, "_", setting2, "_", d,"_cont_qiu.csv"))
}

n_seq <- c(300, 500)
eta_seq <- c(0.5, 1.5)
setting1_seq <- 1
setting2_seq <- 2
d_seq <- c(0.1, 0.2, 0.3, 0.4, 0.5)

for (setting1 in setting1_seq) {
  for (setting2 in setting2_seq) {
    for (n in n_seq) {
      for (eta in eta_seq) {
        for (d in d_seq) {
          sim(
            n = n,
            eta = eta,
            setting1 = setting1,
            setting2 = setting2,
            d = d
          )
        }
      }
    }
  }
}

