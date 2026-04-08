mix.reg.norm.EM <- function(y, X, R, G = 2,
                                 max_iter = 200,
                                 tol = 1e-6,
                                 ridge = 1e-4,
                                 verbose = TRUE){
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  g <- G
  
  #------------------ Stable softmax برداری ------------------#
  softmax_mat <- function(M){
    M <- M - apply(M, 1, max)
    expM <- exp(M)
    expM / rowSums(expM)
  }
  
  #------------------ Normal log-density ------------------#
  logNorm <- function(e, sigma2){
    -0.5*log(2*pi*sigma2) - (e^2)/(2*sigma2)
  }
  
  #------------------ Initialization ------------------#
  km <- kmeans(y, g)
  
  beta <- matrix(0, p, g)
  sigma2 <- rep(var(y), g)
  alpha <- matrix(0, q, g-1)
  
  for(j in 1:g){
    idx <- which(km$cluster == j)
    if(length(idx) > p){
      beta[,j] <- solve(crossprod(X[idx, ]) + ridge*diag(p),
                        crossprod(X[idx, ], y[idx]))
    }
  }
  
  loglik_trace <- numeric(max_iter)
  
  #================ EM =================#
  for(iter in 1:max_iter){
    
    #------------------ E-step  ------------------#
    mu_mat <- X %*% beta     # n x g
    dens_mat <- sapply(1:g, function(j) exp(logNorm(y - mu_mat[,j], sigma2[j])))
    
    pi_mat <- softmax_mat(cbind(R %*% alpha, 0))  # n x g
    tau <- pi_mat * dens_mat
    tau <- pmax(tau, 1e-12)
    tau <- tau / rowSums(tau)
    
    #------------------ M-step ------------------#
    for(j in 1:g){
      w <- tau[,j]
      W <- diag(w)
      beta[,j] <- solve(t(X) %*% W %*% X + ridge*diag(p), t(X) %*% W %*% y)
      
      e_j <- y - X %*% beta[,j]
      sigma2[j] <- sum(w * e_j^2) / sum(w)
      sigma2[j] <- max(sigma2[j], 1e-6)
    }
    
    #------------------ M-step alpha ------------------#
    loglik_alpha <- function(alpha_vec){
      alpha_mat <- matrix(alpha_vec, nrow=q, ncol=g-1)
      pi_mat <- softmax_mat(cbind(R %*% alpha_mat, 0))
      -sum(tau * log(pmax(pi_mat,1e-12)))
    }
    
    opt <- optim(as.vector(alpha), loglik_alpha, method="BFGS")
    alpha <- matrix(opt$par, nrow=q, ncol=g-1)
    
    #------------------ Log-likelihood  ------------------#
    loglik <- sum(log(rowSums(pi_mat * dens_mat)))
    loglik_trace[iter] <- loglik
    
    if(verbose) cat("Iter:", iter, " LogLik:", round(loglik,4), "\n")
    
    if(iter>1 && abs(loglik - loglik_trace[iter-1]) < tol) break
  }
  
  loglik_trace <- loglik_trace[1:iter]
  
  #------------------ Model selection ------------------#
  m <- g*(p+1) + (g-1)*q
  AIC <- -2*loglik + 2*m
  BIC <- -2*loglik + log(n)*m
  EDC = -2 * loglik + 0.2 * sqrt(n) * m
  CAIC = -2 * loglik + 2 * n * m / (n - m - 1)
  ABIC = -2 * loglik + m * log((n + 2) / 24)
  
  #------------------ Output ------------------#
  list(
    beta = beta,
    sigma2 = sigma2,
    alpha = alpha,
    tau = tau,
    loglik = loglik,
    loglik_trace = loglik_trace,
    AIC = AIC,
    BIC = BIC,
    EDC=EDC,
    CAIC=CAIC,
    ABIC=ABIC,
    clusters = apply(tau,1,which.max),
    npar = m
  )
}
