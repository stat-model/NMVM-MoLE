############################################################
# Robust and Complex MoLE Simulation 
############################################################

rm(list = ls())

library(sn)
library(mclust)
library(aricode)
library(MASS)

# ==================== Functions ==================== #
WD.PATH <- paste0(getwd(),"/Functions")

source(paste0(WD.PATH,'/Additional.r'))
source(paste0(WD.PATH,'/GHST.r'))
source(paste0(WD.PATH,'/NMVBS.r'))
source(paste0(WD.PATH,'/NIG.r'))
source(paste0(WD.PATH,'/SMSN_fiting.r'))
source(paste0(WD.PATH,'/normal.r'))

# ----------------- Data Generation ----------------- #
generate_data_complex <- function(n){
  x <- cbind(1, runif(n,-2,2), runif(n,-2,2), runif(n,-2,2), runif(n,-2,2))
  r <- cbind(1, runif(n,-1,3), runif(n,-1,3))
  
  beta <- list(c(-3,2.5,1.5,0,0), c(2,-2,1,0,0), c(0.5,1.5,-2,0,0))
  sigma2 <- c(1.5,2,1.2)
  lambda <- c(2,-1,3)
  
  softmax <- function(v) exp(v)/sum(exp(v))
  
  Z <- y <- numeric(n)
  for(i in 1:n){
    eta <- c(1.5*r[i,2]-2*r[i,3]+0.5*r[i,2]^2,
             -1*r[i,2]+0.5*r[i,3]+0.3*sin(r[i,2]),
             0)
    pi_i <- softmax(eta)
    Z[i] <- sample(1:3,1,prob=pi_i)
    j <- Z[i]
    mu <- sum(x[i,]*beta[[j]])
    y[i] <- rsn(1, xi=mu, omega=sqrt(sigma2[j]), alpha=lambda[j])
  }
  list(y=y, x=x, r=r, z=Z)
}

# ----------------- Metrics ----------------- #
compute_metrics <- function(true_z, est_z, loglik, k, time){
  n <- length(true_z)
  ari <- adjustedRandIndex(true_z, est_z)
  ami <- AMI(true_z, est_z)
  caic <- -2*loglik + k*(log(n)+1)
  c(CAIC=caic, ARI=ari, AMI=ami, TIME=time)
}

# ----------------- Fitting Functions ----------------- #
fit_gaussian <- function(y,x,r,g=3){
  tryCatch({
    t0 <- Sys.time()
    fit <- mix.reg.norm.EM(y=y, X=x, R=r, G=g, verbose=FALSE)
    t1 <- Sys.time()
    list(class=fit$clusters, loglik=fit$loglik,
         npar=g*(ncol(x)+1)+(g-1)*ncol(r),
         time=as.numeric(t1-t0))
  }, error=function(e) NULL)
}

fit_SN <- function(y,x,r,g=3){
  tryCatch({
    t0 <- Sys.time()
    fit <- mix.reg.SMSN.EM(y=y, X=x, R=r, G=g, family="Skew.n", verbose=FALSE)
    t1 <- Sys.time()
    list(class=fit$clusters, loglik=fit$loglik,
         npar=g*(ncol(x)+2)+(g-1)*ncol(r),
         time=as.numeric(t1-t0))
  }, error=function(e) NULL)
}

fit_kmeans <- function(y,x,g=3){
  t0 <- Sys.time()
  km <- kmeans(cbind(y,x), centers=g)
  t1 <- Sys.time()
  list(class=km$cluster, loglik=NA, npar=0, time=as.numeric(t1-t0))
}

fit_CWM <- function(y,x,g=3){
  t0 <- Sys.time()
  fit <- mclust::Mclust(cbind(y,x), G=g)
  t1 <- Sys.time()
  list(class=fit$classification, loglik=fit$loglik, npar=fit$df, time=as.numeric(t1-t0))
}

fit_NIG <- function(y,x,r,g=3){
  t0 <- Sys.time()
  fit <- mix.NIG.MoE.EM(y=y,X=x,R=r,G=g,verbose=FALSE)
  t1 <- Sys.time()
  list(class=fit$clusters, loglik=fit$loglik,
       npar=g*(ncol(x)+3)+(g-1)*ncol(r),
       time=as.numeric(t1-t0))
}

fit_GHST <- function(y,x,r,g=3){
  t0 <- Sys.time()
  fit <- mix.GHST.MoE.EM(y=y,X=x,R=r,G=g,verbose=FALSE)
  t1 <- Sys.time()
  list(class=fit$clusters, loglik=fit$loglik,
       npar=g*(ncol(x)+3)+(g-1)*ncol(r),
       time=as.numeric(t1-t0))
}

fit_NMVBS <- function(y,x,r,g=3){
  t0 <- Sys.time()
  fit <- mix.NMVBS.MoE.EM(y=y,X=x,R=r,G=g,verbose=FALSE)
  t1 <- Sys.time()
  list(class=fit$clusters, loglik=fit$loglik,
       npar=g*(ncol(x)+4)+(g-1)*ncol(r),
       time=as.numeric(t1-t0))
}

# ----------------- Simulation ----------------- #
run_simulation_normal <- function(R=300, n=500){
  methods <- c("kmeans","CWM","Gaussian","SN","NIG","GHST","NMVBS")
  res <- array(NA, c(R, length(methods), 4),
               dimnames=list(NULL, methods, c("CAIC","ARI","AMI","TIME")))
  
  for(rp in 1:R){
    cat("Replication:", rp, "\n")
    data <- generate_data_complex(n)
    fits <- list(
      kmeans = fit_kmeans(data$y,data$x),
      CWM    = fit_CWM(data$y,data$x),
      Gaussian = fit_gaussian(data$y,data$x,data$r),
      SN       = fit_SN(data$y,data$x,data$r),
      NIG      = fit_NIG(data$y,data$x,data$r),
      GHST     = fit_GHST(data$y,data$x,data$r),
      NMVBS    = fit_NMVBS(data$y,data$x,data$r)
    )
    
    for(m in 1:length(methods)){
      f <- fits[[m]]
      if(!is.null(f)){
        res[rp,m,] <- compute_metrics(data$z,f$class,
                                      ifelse(is.na(f$loglik),0,f$loglik),
                                      f$npar,f$time)
      }
    }
  }
  return(res)
}

# ----------------- Summary ----------------- #
summarize_results <- function(res){
  methods <- dimnames(res)[[2]]
  for(m in methods){
    cat("\n====================\nMethod:", m,"\n====================\n")
    means <- colMeans(res[,m,], na.rm=TRUE)
    sds <- apply(res[,m,],2,sd,na.rm=TRUE)
    print(round(cbind(mean=means, sd=sds),3))
  }
}

# ----------------- Run ----------------- #
results <- run_simulation_normal(R=5)
summarize_results(results)
write.csv(results, "simulation3_results.csv")