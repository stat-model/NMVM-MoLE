############################################################
# Robust MoLE Simulation Framework
############################################################

#========================================================
# Clear environment & load packages
#========================================================
rm(list = ls())

library(sn)
library(mclust)
library(aricode)
library(MASS)

#========================================================
# Functions directory
#========================================================
WD.PATH <- paste0(getwd(), "/Functions")

source(file.path(WD.PATH, 'Additional.r'))
source(file.path(WD.PATH, 'GHST.r'))
source(file.path(WD.PATH, 'NMVBS.r'))
source(file.path(WD.PATH, 'NIG.r'))
source(file.path(WD.PATH, 'SMSN_fiting.r'))
source(file.path(WD.PATH, 'normal.r'))
source(file.path(WD.PATH, 'Trimmed.r'))
source(file.path(WD.PATH, 'Contaminated.r'))
source(file.path(WD.PATH, 'RobustBI.r'))

#========================================================
# Utility: Stable softmax
#========================================================
softmax_stable <- function(v){
  # Returns softmax probabilities in a numerically stable way
  exp_v <- exp(v - max(v))
  exp_v / sum(exp_v)
}

#========================================================
# Data generation
#========================================================
generate_data_robust <- function(n, scenario = "S1",
                                 nu = 4, eps = 0.1){
  #' Generate robust mixture-of-experts data
  #' 
  #' @param n sample size
  #' @param scenario "S1" (baseline), "S2" (heavy-tailed), "S3" (leverage contamination)
  #' @param nu degrees of freedom for heavy-tailed (S2)
  #' @param eps contamination proportion (S3)
  
  # Covariates
  x1 <- runif(n, -1, 5)
  r1 <- runif(n, -2, 3)
  
  x <- cbind(1, x1)
  r <- cbind(1, r1)
  
  # True parameters
  beta <- list(c(-4,2), c(3,2), c(0.5,-2))
  sigma2 <- c(2,3,1)
  lambda <- c(3,2,-4)
  alpha <- cbind(c(0.5,-2), c(-1,0.6))
  
  y <- z <- numeric(n)
  
  for(i in 1:n){
    eta <- c(r[i,] %*% alpha, 0)
    pi_i <- softmax_stable(eta)
    
    # component assignment
    z[i] <- sample(1:3,1,prob=pi_i)
    j <- z[i]
    mu <- sum(x[i,]*beta[[j]])
    
    # scenario-specific response
    if(scenario=="S1"){  # skew-normal
      y[i] <- rsn(1, xi=mu, omega=sqrt(sigma2[j]), alpha=lambda[j])
    }
    if(scenario=="S2"){  # skew-t heavy-tailed
      y[i] <- rst(1, xi=mu, omega=sqrt(sigma2[j]), alpha=lambda[j], nu=nu)
    }
    if(scenario=="S3"){  # leverage contamination
      if(runif(1)<eps){
        x[i,2] <- runif(1,10,20)
      }
      mu <- sum(x[i,]*beta[[j]])
      y[i] <- rsn(1, xi=mu, omega=sqrt(sigma2[j]), alpha=lambda[j])
    }
  }
  
  return(list(y=y, x=x, r=r, z=z))
}

#========================================================
# Metrics computation
#========================================================
compute_metrics <- function(true_z, est_z, loglik, npar){
  ari <- adjustedRandIndex(true_z, est_z)
  ami <- AMI(true_z, est_z)
  n <- length(true_z)
  caic <- -2*loglik + npar*(log(n)+1)
  c(CAIC=caic, ARI=ari, AMI=ami)
}

#========================================================
# Run simulation
#========================================================
run_robust_simulation <- function(R = 200, n = 500){
  # Run Monte Carlo simulation for robust MoLE models
  scenarios <- list(
    S1   = list(type="S1"),
    S2_6 = list(type="S2", nu=6),
    S2_4 = list(type="S2", nu=4),
    S2_2 = list(type="S2", nu=2),
    S3_5  = list(type="S3", eps=0.05),
    S3_10 = list(type="S3", eps=0.10),
    S3_15 = list(type="S3", eps=0.15)
  )
  
  methods <- c("Gaussian","Trimmed","Contaminated",
               "RobustBI","Skew","NIG","GHST","NMVBS")
  
  results <- list()
  
  for(sc in names(scenarios)){
    cat("Scenario:", sc, "\n")
    res_mat <- array(NA, dim = c(R, length(methods), 3),
                     dimnames = list(NULL, methods, c("CAIC","ARI","AMI")))
    
    for(rp in 1:R){
      sc_info <- scenarios[[sc]]
      
      data <- switch(sc_info$type,
                     S1 = generate_data_robust(n, "S1"),
                     S2 = generate_data_robust(n, "S2", nu=sc_info$nu),
                     S3 = generate_data_robust(n, "S3", eps=sc_info$eps))
      
      y <- data$y; x <- data$x; r <- data$r; z <- data$z
      
      # fit models
      fits <- list(
        Gaussian = tryCatch(mix.reg.norm.EM(y, x, r, G = 3, verbose = FALSE), error=function(e) NULL),
        Trimmed = tryCatch(fit_trimmed_mole(y,x,r,g=3), error=function(e) NULL),
        Contaminated = tryCatch(fit_contaminated_mole(y,x,r,G=3), error=function(e) NULL),
        RobustBI = tryCatch(fit_robust_bi_mole(y,x,r,G=3), error=function(e) NULL),
        Skew = tryCatch(mix.reg.SMSN.EM(y, x, r, G = 3, family="Skew.n", verbose=FALSE), error=function(e) NULL),
        NIG = tryCatch(mix.NIG.MoE.EM(y, x, r, G = 3, verbose=FALSE), error=function(e) NULL),
        GHST = tryCatch(mix.GHST.MoE.EM(y,x,r, G = 3, verbose=FALSE), error=function(e) NULL),
        NMVBS = tryCatch(mix.NMVBS.MoE.EM(y, x, r, G = 3, verbose=FALSE), error=function(e) NULL)
      )
      
      # compute metrics
      for(m in methods){
        f <- fits[[m]]
        if(!is.null(f)){
          res_mat[rp, m, ] <- compute_metrics(z, f$clusters, f$loglik, f$npar)
        }
      }
    }
    results[[sc]] <- res_mat
  }
  return(results)
}

#========================================================
# Summarize results
#========================================================
summarize_results <- function(results){
  #' Print mean ± SD for CAIC, ARI, AMI
  for(sc in names(results)){
    cat("\n====================\nScenario:", sc, "\n====================\n")
    arr <- results[[sc]]
    for(m in dimnames(arr)[[2]]){
      vals <- arr[, m, ]
      means <- colMeans(vals, na.rm=TRUE)
      sds <- apply(vals, 2, sd, na.rm=TRUE)
      cat("\n", m, ":\n")
      cat(sprintf("CAIC: %.3f (%.3f)\n", means["CAIC"], sds["CAIC"]))
      cat(sprintf("ARI : %.3f (%.3f)\n", means["ARI"], sds["ARI"]))
      cat(sprintf("AMI : %.3f (%.3f)\n", means["AMI"], sds["AMI"]))
    }
  }
}

#========================================================
# Run simulation & summarize
#========================================================
results <- run_robust_simulation(R)
summarize_results(results)
write.csv(results, "simulation2_results.csv")