mix.NMVBS.MoE.EM <- function(y, X, R, G = 2,
                                  max_iter = 200,
                                  tol = 1e-6,
                                  ridge = 1e-4,
                                  verbose = TRUE) {
  library(ghyp)
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  
  #------------------ Softmax برداری ------------------#
  softmax_mat <- function(M) {
    M <- M - apply(M, 1, max)
    expM <- exp(M)
    expM / rowSums(expM)
  }
  
  #------------------ GH log-density ------------------#
  logGH_vec <- function(e, sigma2, lambda, kappa, alpha2) {
    chi <- alpha2^(-2)
    psi <- alpha2^(-2)
    obj <- ghyp(lambda = kappa, chi = chi, psi = psi,
                sigma = sqrt(sigma2), gamma = lambda)
    dghyp(e, object = obj, logvalue = TRUE)
  }
  
  #------------------ NMVBS log-density ------------------#
  logNMVBS_vec <- function(e, sigma2, lambda, alpha2) {
    l1 <- logGH_vec(e, sigma2, lambda, 0.5, alpha2)
    l2 <- logGH_vec(e, sigma2, lambda, -0.5, alpha2)
    m <- pmax(l1, l2)
    m + log(0.5 * exp(l1 - m) + 0.5 * exp(l2 - m))
  }
  
  #------------------ Initialization ------------------#
  km <- kmeans(y, G)
  
  beta <- matrix(0, p, G)
  sigma2 <- rep(var(y), G)
  lambda <- rep(0, G)
  alpha2 <- rep(1, G)
  alpha <- matrix(0, q, G - 1)
  
  for (j in 1:G) {
    idx <- which(km$cluster == j)
    if (length(idx) > p) {
      beta[, j] <- solve(crossprod(X[idx, ]) + ridge * diag(p),
                         crossprod(X[idx, ], y[idx]))
    }
  }
  
  loglik_trace <- numeric(max_iter)
  
  #================ EM =================#
  for (iter in 1:max_iter) {
    
    #------------------ E-step ------------------#
    mu_mat <- X %*% beta  # n x G
    dens_mat <- matrix(0, n, G)
    for (j in 1:G) {
      dens_mat[, j] <- exp(logNMVBS_vec(y - mu_mat[, j],
                                        sigma2[j],
                                        lambda[j],
                                        alpha2[j]))
    }
    
    eta_mat <- cbind(R %*% alpha, 0)
    pi_mat <- softmax_mat(eta_mat)
    
    tau_ij <- pi_mat * dens_mat
    tau_ij <- pmax(tau_ij, 1e-12)
    tau_ij <- tau_ij / rowSums(tau_ij)
    
    #------------------ M-step (beta, sigma2, lambda) ------------------#
    for (j in 1:G) {
      W <- diag(tau_ij[, j])
      beta[, j] <- solve(t(X) %*% W %*% X + ridge * diag(p),
                         t(X) %*% W %*% y)
      e_j <- y - X %*% beta[, j]
      w <- tau_ij[, j]
      lambda[j] <- sum(w * e_j) / sum(w)
      lambda[j] <- max(min(lambda[j], 10), -10)
      sigma2[j] <- sum(w * (e_j - lambda[j])^2) / sum(w)
      sigma2[j] <- min(max(sigma2[j], 1e-4), 1e4)
    }
    
    #------------------ M-step (alpha2) ------------------#
    for (j in 1:G) {
      obj_a2 <- function(a2) {
        if (a2 <= 0) return(1e10)
        sum_tau_log <- sum(tau_ij[, j] *
                             logNMVBS_vec(y - X %*% beta[, j],
                                          sigma2[j],
                                          lambda[j],
                                          a2))
        -sum_tau_log
      }
      opt <- optim(alpha2[j], obj_a2,
                   method = "L-BFGS-B",
                   lower = 1e-3, upper = 50)
      alpha2[j] <- opt$par
    }
    
    #------------------ M-step (alpha gating) ------------------#
    loglik_alpha <- function(alpha_vec) {
      alpha_mat <- matrix(alpha_vec, nrow = q, ncol = G - 1)
      eta_mat <- cbind(R %*% alpha_mat, 0)
      pi_mat <- softmax_mat(eta_mat)
      -sum(tau_ij * log(pmax(pi_mat, 1e-12)))
    }
    opt <- optim(as.vector(alpha), loglik_alpha,
                 method = "BFGS", control = list(maxit = 100))
    alpha <- matrix(opt$par, nrow = q, ncol = G - 1)
    
    #------------------ Log-likelihood ------------------#
    loglik <- sum(log(rowSums(pi_mat * dens_mat)))
    loglik_trace[iter] <- loglik
    if (verbose) cat("Iter:", iter, "LogLik:", loglik, "\n")
    if (iter > 1 && abs(loglik - loglik_trace[iter-1]) < tol) break
  }
  
  loglik_trace <- loglik_trace[1:iter]
  m <- G * (p + 3 + 1) + (G - 1) * q
  
  list(
    beta = beta,
    sigma2 = sigma2,
    lambda = lambda,
    nu = alpha2,
    alpha = alpha,
    tau = tau_ij,
    loglik = loglik,
    loglik_trace = loglik_trace,
    AIC = -2*loglik + 2*m,
    BIC = -2*loglik + log(n)*m,
    EDC = -2*loglik + 0.2*sqrt(n)*m,
    CAIC = -2*loglik + 2*n*m/(n-m-1),
    ABIC = -2*loglik + m*log((n+2)/24),
    clusters = apply(tau_ij, 1, which.max),
    npar = m
  )
}

