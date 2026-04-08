mix.GHST.MoE.EM <- function(y, X, R, G = 2,
                                 max_iter = 200,
                                 tol = 1e-6,
                                 ridge = 1e-4,
                                 verbose = TRUE) {
  library(ghyp)
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  
  #------------------ Stable softmax (vectorized) ------------------#
  softmax_mat <- function(M) {
    M <- M - apply(M, 1, max)
    expM <- exp(M)
    expM / rowSums(expM)
  }
  
  #------------------ Fast GHST density ------------------#
  dGHST_vec <- function(y, mu, sigma, gamma, nu) {
    obj <- ghyp(lambda = -nu/2, chi = nu, psi = 0,
                sigma = sigma, gamma = gamma, mu = 0)
    dghyp(y - mu, object = obj, logvalue = FALSE)
  }
  
  #------------------ Initialization ------------------#
  km <- kmeans(y, G)
  beta <- matrix(0, p, G)
  sigma <- rep(sd(y), G)
  gamma <- rep(0, G)
  nu <- rep(10, G)
  alpha <- matrix(0, q, G - 1)
  
  for (j in 1:G) {
    idx <- which(km$cluster == j)
    if (length(idx) > p) {
      XtX <- t(X[idx, ]) %*% X[idx, ]
      beta[, j] <- solve(XtX + ridge * diag(p)) %*% t(X[idx, ]) %*% y[idx]
    }
  }
  
  loglik_trace <- numeric(max_iter)
  
  #------------------ EM Algorithm ------------------#
  for (iter in 1:max_iter) {
    
    # ================= E-step =================
    
    # mu for all (n × G)
    mu_mat <- X %*% beta
    
    # densities
    dens_mat <- matrix(0, n, G)
    for (j in 1:G) {
      dens_mat[, j] <- dGHST_vec(y, mu_mat[, j], sigma[j], gamma[j], nu[j])
    }
    
    # gating probabilities
    eta_mat <- cbind(R %*% alpha, 0)
    pi_mat <- softmax_mat(eta_mat)
    
    # responsibilities
    tau_ij <- pi_mat * dens_mat
    tau_ij <- pmax(tau_ij, 1e-12)
    tau_ij <- tau_ij / rowSums(tau_ij)
    
    # ================= M-step =================
    
    # Update beta
    for (j in 1:G) {
      W <- diag(tau_ij[, j])
      XtW <- t(X) %*% W
      beta[, j] <- solve(XtW %*% X + ridge * diag(p)) %*% XtW %*% y
    }
    
    # Update sigma and gamma
    for (j in 1:G) {
      mu_j <- mu_mat[, j]
      res_j <- y - mu_j
      w <- tau_ij[, j]
      
      gamma[j] <- sum(w * res_j) / sum(w)
      gamma[j] <- max(min(gamma[j], 3*sigma[j]), -3*sigma[j])
      
      sigma[j] <- sqrt(sum(w * (res_j - gamma[j])^2) / sum(w))
      sigma[j] <- max(sigma[j], 1e-4)
    }
    
    # Update nu
    for (j in 1:G) {
      obj_nu <- function(nu_j) {
        if (nu_j <= 0) return(1e10)
        mu_j <- mu_mat[, j]
        dens <- dGHST_vec(y, mu_j, sigma[j], gamma[j], nu_j)
        -sum(tau_ij[, j] * log(pmax(dens, 1e-12)))
      }
      opt <- optim(nu[j], obj_nu, method = "L-BFGS-B", lower = 0.5, upper = 50)
      nu[j] <- opt$par
    }
    
    # Update alpha
    loglik_alpha <- function(alpha_vec) {
      alpha_mat <- matrix(alpha_vec, nrow = q, ncol = G - 1)
      eta_mat <- cbind(R %*% alpha_mat, 0)
      pi_mat <- softmax_mat(eta_mat)
      -sum(tau_ij * log(pmax(pi_mat, 1e-12)))
    }
    
    opt <- optim(as.vector(alpha), loglik_alpha,
                 method = "BFGS", control = list(maxit = 100))
    alpha <- matrix(opt$par, nrow = q, ncol = G - 1)
    
    # ================= Log-likelihood =================
    
    loglik <- sum(log(rowSums(pi_mat * dens_mat)))
    loglik_trace[iter] <- loglik
    
    if (verbose) cat("Iter:", iter, "LogLik:", loglik, "\n")
    
    if (iter > 1 && abs(loglik - loglik_trace[iter-1]) < tol) break
  }
  
  loglik_trace <- loglik_trace[1:iter]
  
  m <- G * (p + 3) + (G - 1) * q
  
  list(
    beta = beta,
    sigma2 = sigma,
    lambda = gamma,
    nu = nu,
    alpha = alpha,
    tau = tau_ij,
    loglik = loglik,
    AIC = -2*loglik + 2*m,
    BIC = -2*loglik + log(n)*m,
    EDC = -2 * loglik + 0.2 * sqrt(n) * m,
    CAIC = -2 * loglik + 2 * n * m / (n - m - 1),
    ABIC = -2 * loglik + m * log((n + 2) / 24),
    loglik_trace = loglik_trace,
    clusters = apply(tau_ij, 1, which.max),
    npar = m
  )
}
