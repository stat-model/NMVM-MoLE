############################################################
# Robust MoLE Simulation Framework
# S1: baseline
# S2: heavy-tailed noise
# S3: bad leverage contamination
############################################################

rm(list = ls())

library(sn)
library(mclust)
library(aricode)
library(MASS)

set.seed(12345)

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

if (!dir.exists("results")) dir.create("results")

softmax_stable <- function(v) {
  exp_v <- exp(v - max(v))
  exp_v / sum(exp_v)
}

generate_data_robust <- function(n, scenario = "S1", nu = 4, eps = 0.10, delta = 15) {
  
  x1 <- runif(n, -1, 5)
  r1 <- runif(n, -2, 3)
  
  x <- cbind(1, x1)
  r <- cbind(1, r1)
  
  beta <- list(c(-4, 2), c(3, 2), c(0.5, -2))
  sigma2 <- c(2, 3, 1)
  lambda <- c(3, 2, -4)
  
  alpha <- cbind(c(0.5, -2), c(-1, 0.6))
  
  y <- numeric(n)
  z <- numeric(n)
  
  for (i in 1:n) {
    
    eta <- c(as.numeric(r[i, ] %*% alpha), 0)
    pi_i <- softmax_stable(eta)
    
    z[i] <- sample(1:3, 1, prob = pi_i)
    j <- z[i]
    
    mu <- sum(x[i, ] * beta[[j]])
    
    if (scenario == "S1") {
      
      y[i] <- rsn(
        n = 1,
        xi = mu,
        omega = sqrt(sigma2[j]),
        alpha = lambda[j]
      )
      
    } else if (scenario == "S2") {
      
      y[i] <- rst(
        n = 1,
        xi = mu,
        omega = sqrt(sigma2[j]),
        alpha = lambda[j],
        nu = nu
      )
      
    } else if (scenario == "S3") {
      
      if (runif(1) < eps) {
        
        x[i, 2] <- runif(1, 10, 20)
        
        k <- switch(
          as.character(j),
          "1" = 2,
          "2" = 3,
          "3" = 1
        )
        
        mu_bad <- sum(x[i, ] * beta[[k]]) + delta
        
        y[i] <- rnorm(
          n = 1,
          mean = mu_bad,
          sd = sqrt(sigma2[k])
        )
        
      } else {
        
        y[i] <- rsn(
          n = 1,
          xi = mu,
          omega = sqrt(sigma2[j]),
          alpha = lambda[j]
        )
      }
      
    } else {
      stop("Unknown scenario.")
    }
  }
  
  list(y = y, x = x, r = r, z = z)
}

compute_metrics <- function(true_z, est_z, loglik, npar) {
  
  n <- length(true_z)
  
  ari <- adjustedRandIndex(true_z, est_z)
  ami <- AMI(true_z, est_z)
  caic <- -2 * loglik + npar * (log(n) + 1)
  
  c(CAIC = caic, ARI = ari, AMI = ami)
}

extract_fit_info <- function(f) {
  
  if (is.null(f)) return(NULL)
  
  clusters <- NULL
  if (!is.null(f$clusters)) clusters <- f$clusters
  if (is.null(clusters) && !is.null(f$class)) clusters <- f$class
  if (is.null(clusters) && !is.null(f$classification)) clusters <- f$classification
  
  loglik <- NULL
  if (!is.null(f$loglik)) loglik <- f$loglik
  if (is.null(loglik) && !is.null(f$logLik)) loglik <- f$logLik
  if (is.null(loglik) && !is.null(f$ll)) loglik <- f$ll
  
  npar <- NULL
  if (!is.null(f$npar)) npar <- f$npar
  if (is.null(npar) && !is.null(f$df)) npar <- f$df
  
  if (is.null(clusters) || is.null(loglik) || is.null(npar)) return(NULL)
  
  list(
    clusters = clusters,
    loglik = as.numeric(loglik)[1],
    npar = as.numeric(npar)[1]
  )
}

run_robust_simulation <- function(R = 200, n = 500) {
  
  scenarios <- list(
    S1     = list(type = "S1"),
    S2_6   = list(type = "S2", nu = 6),
    S2_4   = list(type = "S2", nu = 4),
    S2_2   = list(type = "S2", nu = 2),
    S3_2.5 = list(type = "S3", eps = 0.025),
    S3_5   = list(type = "S3", eps = 0.05),
    S3_10  = list(type = "S3", eps = 0.10)
  )
  
  methods <- c(
    "Gaussian",
    "Trimmed",
    "Contaminated",
    "RobustBI",
    "Skew",
    "NIG",
    "GHST",
    "NMVBS"
  )
  
  results <- list()
  
  for (sc in names(scenarios)) {
    
    cat("\nScenario:", sc, "\n")
    
    res_mat <- array(
      NA_real_,
      dim = c(R, length(methods), 3),
      dimnames = list(NULL, methods, c("CAIC", "ARI", "AMI"))
    )
    
    for (rp in 1:R) {
      
      sc_info <- scenarios[[sc]]
      
      data <- switch(
        sc_info$type,
        S1 = generate_data_robust(n, scenario = "S1"),
        S2 = generate_data_robust(n, scenario = "S2", nu = sc_info$nu),
        S3 = generate_data_robust(n, scenario = "S3", eps = sc_info$eps, delta = 15)
      )
      
      y <- data$y
      x <- data$x
      r <- data$r
      z <- data$z
      
      fits <- list(
        Gaussian = tryCatch(
          mix.reg.norm.EM(y = y, X = x, R = r, G = 3, verbose = FALSE),
          error = function(e) NULL
        ),
        Trimmed = tryCatch(
          fit_trimmed_mole(y, x, r, g = 3),
          error = function(e) NULL
        ),
        Contaminated = tryCatch(
          fit_contaminated_mole(y, x, r, G = 3),
          error = function(e) NULL
        ),
        RobustBI = tryCatch(
          fit_robust_bi_mole(y, x, r, G = 3),
          error = function(e) NULL
        ),
        Skew = tryCatch(
          mix.reg.SMSN.EM(y = y, X = x, R = r, G = 3,
                          family = "Skew.n", verbose = FALSE),
          error = function(e) NULL
        ),
        NIG = tryCatch(
          mix.NIG.MoE.EM(y = y, X = x, R = r, G = 3, verbose = FALSE),
          error = function(e) NULL
        ),
        GHST = tryCatch(
          mix.GHST.MoE.EM(y = y, X = x, R = r, G = 3, verbose = FALSE),
          error = function(e) NULL
        ),
        NMVBS = tryCatch(
          mix.NMVBS.MoE.EM(y = y, X = x, R = r, G = 3, verbose = FALSE),
          error = function(e) NULL
        )
      )
      
      for (m in methods) {
        
        info <- extract_fit_info(fits[[m]])
        
        if (!is.null(info)) {
          res_mat[rp, m, ] <- compute_metrics(
            true_z = z,
            est_z = info$clusters,
            loglik = info$loglik,
            npar = info$npar
          )
        }
      }
      
      if (rp %% 10 == 0) {
        cat("Replication:", rp, "of", R, "\n")
      }
    }
    
    results[[sc]] <- res_mat
    
    saveRDS(results, file = "results/temp_results.rds")
    cat("Intermediate results saved for scenario:", sc, "\n")
  }
  
  results
}

summarize_results <- function(results) {
  
  summary_list <- list()
  
  for (sc in names(results)) {
    
    arr <- results[[sc]]
    
    for (m in dimnames(arr)[[2]]) {
      
      vals <- arr[, m, ]
      
      means <- colMeans(vals, na.rm = TRUE)
      sds <- apply(vals, 2, sd, na.rm = TRUE)
      success <- sum(!is.na(vals[, "ARI"]))
      
      summary_list[[length(summary_list) + 1]] <- data.frame(
        Scenario = sc,
        Method = m,
        CAIC_mean = means["CAIC"],
        CAIC_sd = sds["CAIC"],
        ARI_mean = means["ARI"],
        ARI_sd = sds["ARI"],
        AMI_mean = means["AMI"],
        AMI_sd = sds["AMI"],
        Successful_fits = success
      )
    }
  }
  
  do.call(rbind, summary_list)
}

# Test first if needed:
# results <- run_robust_simulation(R = 2, n = 100)

results <- run_robust_simulation(R = 200, n = 500)

summary_results <- summarize_results(results)

print(summary_results)

saveRDS(results, "results/simulation_4_2_full_results.rds")

write.csv(
  summary_results,
  "results/simulation_4_2_summary.csv",
  row.names = FALSE
)
