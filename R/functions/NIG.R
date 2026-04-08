mix.NIG.MoE.EM <- function(y, X, R, G = 2,
                                max_iter = 200,
                                tol = 1e-6,
                                ridge = 1e-4,
                                verbose = TRUE) {
  
  library(ghyp)
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  
  #------------------ Stable softmax ------------------#
  softmax_mat <- function(eta_mat) {
    eta_mat <- eta_mat - apply(eta_mat, 1, max)
    exp_eta <- exp(eta_mat)
    exp_eta / rowSums(exp_eta)
  }
  
  #------------------ Initialization ------------------#
  km <- kmeans(y, G)
  
  beta <- matrix(0, p, G)
  sigma <- rep(sqrt(var(y)), G)
  gamma <- rep(0, G)
  nu <- rep(5, G)
  alpha <- matrix(0, q, G - 1)
  
  for (j in 1:G) {
    idx <- which(km$cluster == j)
    if (length(idx) > p) {
      XtX <- t(X[idx, ]) %*% X[idx, ]
      beta[, j] <- solve(XtX + ridge * diag(p)) %*% t(X[idx, ]) %*% y[idx]
    }
  }
  
  loglik_trace <- c()
  
  #=============================#
  # EM Algorithm
  #=============================#
  for (iter in 1:max_iter) {
    
    #-------------------------#
    # E-step (vectorized)
    #-------------------------#
    
    # gating
    eta <- cbind(R %*% alpha, 0)
    pi_mat <- softmax_mat(eta)
    
    dens_mat <- matrix(0, n, G)
    
    for (j in 1:G) {
      mu_j <- X %*% beta[, j]
      e_j <- y - mu_j
      
      obj <- NIG(mu = 0, sigma = sigma[j], gamma = gamma[j], psi = nu[j], chi = 1)
      
      logdens <- dghyp(e_j, object = obj, logvalue = TRUE)
      dens_mat[, j] <- exp(pmax(logdens, -700))
    }
    
    tau_ij <- pi_mat * dens_mat
    tau_ij <- pmax(tau_ij, 1e-12)
    tau_ij <- tau_ij / rowSums(tau_ij)
    
    #-------------------------#
    # M-step: beta
    #-------------------------#
    for (j in 1:G) {
      W <- diag(tau_ij[, j])
      XtW <- t(X) %*% W
      beta[, j] <- solve(XtW %*% X + ridge * diag(p)) %*% XtW %*% y
    }
    
    #-------------------------#
    # M-step: NIG params
    # (فقط vectorize داخل obj)
    #-------------------------#
    for (j in 1:G) {
      
      obj_fun <- function(par) {
        sigma_j <- par[1]
        gamma_j <- par[2]
        nu_j <- par[3]
        
        if (sigma_j <= 0 || nu_j <= 0) return(1e10)
        
        mu_j <- X %*% beta[, j]
        e_j <- y - mu_j
        
        obj <- NIG(mu = 0, sigma = sigma_j, gamma = gamma_j, psi = nu_j, chi = 1)
        logdens <- dghyp(e_j, object = obj, logvalue = TRUE)
        
        ll <- sum(tau_ij[, j] * pmax(logdens, -700))
        return(-ll)
      }
      
      opt <- optim(c(sigma[j], gamma[j], nu[j]), obj_fun,
                   method = "L-BFGS-B",
                   lower = c(1e-4, -10, 0.1),
                   upper = c(1e4, 10, 50))
      
      sigma[j] <- opt$par[1]
      gamma[j] <- opt$par[2]
      nu[j] <- opt$par[3]
    }
    
    #-------------------------#
    # M-step: alpha
    #-------------------------#
    loglik_alpha <- function(alpha_vec) {
      alpha_mat <- matrix(alpha_vec, nrow = q, ncol = G - 1)
      eta <- cbind(R %*% alpha_mat, 0)
      pi_mat <- softmax_mat(eta)
      -sum(tau_ij * log(pmax(pi_mat, 1e-12)))
    }
    
    opt <- optim(as.vector(alpha), loglik_alpha,
                 method = "BFGS",
                 control = list(maxit = 100))
    
    alpha <- matrix(opt$par, nrow = q, ncol = G - 1)
    
    #-------------------------#
    # Log-likelihood
    #-------------------------#
    loglik <- sum(log(rowSums(pi_mat * dens_mat)))
    loglik_trace <- c(loglik_trace, loglik)
    
    if (verbose) cat("Iter:", iter, "LogLik:", loglik, "\n")
    
    if (iter > 1 && abs(loglik - loglik_trace[iter - 1]) < tol) break
  }
  
  m <- G * (p + 3) + (G - 1) * q
  
  list(
    beta = beta,
    sigma = sigma,
    gamma = gamma,
    nu = nu,
    alpha = alpha,
    tau = tau_ij,
    loglik = loglik,
    AIC = -2 * loglik + 2 * m,
    BIC = -2 * loglik + log(n) * m,
    EDC = -2 * loglik + 0.2 * sqrt(n) * m,
    CAIC = -2 * loglik + 2 * n * m / (n - m - 1),
    ABIC = -2 * loglik + m * log((n + 2) / 24),
    loglik_trace = loglik_trace,
    clusters = apply(tau_ij, 1, which.max),
    npar = m
  )
}

