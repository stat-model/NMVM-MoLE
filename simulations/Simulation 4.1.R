# =========================================
# NMVM-MoLE Parameter Recovery Simulation
# =========================================

rm(list = ls())

# ------------------- Libraries ------------------- #
library(sn)
library(ghyp)

# ------------------- Source Functions ------------------- #
WD.PATH <- paste0(getwd(), "/Functions")
source(paste0(WD.PATH, '/Additional.r'))
source(paste0(WD.PATH, '/GHST.r'))
source(paste0(WD.PATH, '/NMVBS.r'))

# Create folder to store results
if(!dir.exists("results")) dir.create("results")

# ------------------- True Parameters ------------------- #
true_param <- list(
  beta1 = c(-1, -2, -3),
  beta2 = c(3, 7, 2),
  sigma2 = c(1, 2),
  lambda = c(3, -2),
  nu = c(3, 5),
  tau = c(-1, 0.5, 3)
)

true_vec <- c(true_param$beta1,
              true_param$beta2,
              true_param$sigma2,
              true_param$lambda,
              true_param$nu,
              true_param$tau)

# ------------------- Data Generation ------------------- #
generate_data <- function(n, family = "GHST"){
  X <- cbind(1, runif(n, -1, 4), runif(n, -3, 3))
  R <- cbind(1, runif(n, -2, 1), runif(n, -1, 1))
  beta <- cbind(true_param$beta1, true_param$beta2)
  
  theta <- if(family == "GHST") matrix(true_param$nu, nrow = 1) else matrix(c(1.5, 0.8), nrow = 1)
  
  sim <- r.mix.NMVM(
    n = n,
    X = X,
    beta = beta,
    sigma2 = true_param$sigma2,
    lambda = true_param$lambda,
    theta = theta,
    alpha = matrix(true_param$tau, ncol = 1),
    R = R,
    family = family
  )
  list(y = sim$y, x = X, r = R, z = sim$class)
}

# ------------------- Reorder Components ------------------- #
reorder_components <- function(fit, true_beta){
  # Reorder mixture components based on proximity to true beta
  d1 <- sum((fit$beta[,1] - true_beta[,1])^2)
  d2 <- sum((fit$beta[,1] - true_beta[,2])^2)
  if(d2 < d1){
    fit$beta   <- fit$beta[,c(2,1)]
    fit$sigma2 <- fit$sigma2[c(2,1)]
    fit$lambda <- fit$lambda[c(2,1)]
    fit$nu     <- fit$nu[c(2,1)]
  }
  fit
}

# ------------------- Extract Parameters ------------------- #
extract_param <- function(fit){
  c(fit$beta[,1], fit$beta[,2],
    fit$sigma2,
    fit$lambda,
    fit$nu,
    as.vector(fit$alpha))
}

# ------------------- Compute Metrics ------------------- #
compute_metrics <- function(est_mat, true_vec){
  # Bias, standard deviation, RMSE, absolute relative bias
  Bias = colMeans(est_mat - matrix(true_vec, nrow=nrow(est_mat), ncol=length(true_vec), byrow=TRUE))
  STD  = apply(est_mat, 2, sd)
  RMSE = sqrt(colMeans((est_mat - matrix(true_vec, nrow=nrow(est_mat), ncol=length(true_vec), byrow=TRUE))^2))
  ARB  = colMeans(abs((est_mat - matrix(true_vec, nrow=nrow(est_mat), ncol=length(true_vec), byrow=TRUE)) / 
                        pmax(abs(matrix(true_vec, nrow=nrow(est_mat), ncol=length(true_vec), byrow=TRUE)),1e-6)))
  data.frame(Bias=Bias, STD=STD, RMSE=RMSE, ARB=ARB)
}

# ------------------- Simulation ------------------- #
run_simulation <- function(family = "GHST", R = 500, n = 200){
  est_mat <- matrix(NA, R, length(true_vec))
  success <- 0
  for(rep in 1:R){
    data <- generate_data(n, family)
    true_beta <- cbind(true_param$beta1, true_param$beta2)
    
    fit <- tryCatch({
      if(family == "GHST"){
        mix.GHST.MoE.EM(y = data$y, X = data$x, R = data$r, G = 2, tol = 1e-6, verbose = FALSE)
      } else {
        mix.NMVBS.MoE.EM(y = data$y, X = data$x, R = data$r, G = 2, tol = 1e-6, verbose = FALSE)
      }
    }, error = function(e) NULL)
    
    if(!is.null(fit)){
      fit <- reorder_components(fit, true_beta)
      success <- success + 1
      est_mat[success,] <- extract_param(fit)
    }
    
    if(rep %% 50 == 0) cat("Replication:", rep, "\n")
  }
  
  est_mat <- est_mat[1:success, , drop=FALSE]
  metrics <- compute_metrics(est_mat, true_vec)
  
  # Save results to CSV
  write.csv(metrics, paste0("results/metrics_", family, "_n", n, ".csv"), row.names = TRUE)
  metrics
}

# ------------------- Run  Simulation ------------------- #
results_GHST  <- run_simulation("GHST", R=500, n=200)
results_NMVBS <- run_simulation("NMVBS", R=500, n=200)

cat("=== GHST-MoLE Metrics ===\n"); print(results_GHST)
cat("\n=== NMVBS-MoLE Metrics ===\n"); print(results_NMVBS)
