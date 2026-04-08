mix.VG.MoE.EM <- function(y, X, R, G = 2,
                               max_iter = 200, 
                               tol = 1e-6,
                               ridge = 1e-4,
                               verbose = TRUE) {
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  
  library(ghyp)
  
  #------------------ Stable softmax  ------------------#
  softmax_mat <- function(M){
    M <- M - apply(M, 1, max)
    expM <- exp(M)
    expM / rowSums(expM)
  }
  
  #------------------ VG log-density ------------------#
  logVG_vec <- function(e, sigma2, lambda, kappa, psi){
    obj <- VG(lambda = kappa,
              psi = psi,
              mu = 0,
              sigma = sqrt(sigma2),
              gamma = lambda)
    dghyp(e, object = obj, logvalue = TRUE)
  }
  
  #------------------ Initialization ------------------#
  km <- kmeans(y, G)
  
  beta <- matrix(0, p, G)
  sigma2 <- rep(var(y), G)
  lambda <- rep(0, G)
  nu <- matrix(1, 2, G)  # kappa, psi
  alpha <- matrix(0, q, G-1)
  
  for(j in 1:G){
    idx <- which(km$cluster == j)
    if(length(idx) > p){
      beta[, j] <- solve(crossprod(X[idx, ]) + ridge*diag(p),
                         crossprod(X[idx, ], y[idx]))
    }
  }
  
  loglik_trace <- numeric(max_iter)
  
  #================ EM =================#
  for(iter in 1:max_iter){
    
    #------------------ E-step (برداری) ------------------#
    mu_mat <- X %*% beta     # n x G
    dens_mat <- matrix(0, n, G)
    
    for(j in 1:G){
      dens_mat[, j] <- logVG_vec(y - mu_mat[, j], sigma2[j], lambda[j],
                                 nu[1, j], nu[2, j]) |> exp()
    }
    
    pi_mat <- softmax_mat(cbind(R %*% alpha, 0))
    
    tau_ij <- pi_mat * dens_mat
    tau_ij <- pmax(tau_ij, 1e-12)
    tau_ij <- tau_ij / rowSums(tau_ij)
    
    #------------------ M-step (beta, sigma2, lambda, nu) ------------------#
    for(j in 1:G){
      w <- tau_ij[, j]
      W <- diag(w)
      
      beta[, j] <- solve(t(X) %*% W %*% X + ridge*diag(p),
                         t(X) %*% W %*% y)
      
      e_j <- y - X %*% beta[, j]
      
      lambda[j] <- sum(w * e_j)/sum(w)
      lambda[j] <- max(min(lambda[j], 10), -10)
      
      sigma2[j] <- sum(w * (e_j - lambda[j])^2)/sum(w)
      sigma2[j] <- min(max(sigma2[j], 1e-4), 1e4)
      
      # VG parameters
      obj_fun <- function(par){
        kappa_j <- par[1]
        psi_j <- par[2]
        if(kappa_j <= 0 || psi_j <= 0) return(1e10)
        -sum(w * logVG_vec(e_j, sigma2[j], lambda[j], kappa_j, psi_j))
      }
      
      opt <- optim(c(nu[1, j], nu[2, j]), obj_fun,
                   method="L-BFGS-B", lower=c(1e-4, 1e-4), upper=c(50, 50))
      nu[, j] <- opt$par
    }
    
    #------------------ M-step (alpha gating) ------------------#
    loglik_alpha <- function(alpha_vec){
      alpha_mat <- matrix(alpha_vec, nrow=q, ncol=G-1)
      pi_mat <- softmax_mat(cbind(R %*% alpha_mat, 0))
      -sum(tau_ij * log(pmax(pi_mat, 1e-12)))
    }
    
    opt <- optim(as.vector(alpha), loglik_alpha,
                 method="BFGS", control=list(maxit=100))
    alpha <- matrix(opt$par, nrow=q, ncol=G-1)
    
    #------------------ Log-likelihood (برداری) ------------------#
    loglik <- sum(log(rowSums(pi_mat * dens_mat)))
    loglik_trace[iter] <- loglik
    
    if(verbose) cat("Iter:", iter, " LogLik:", loglik, "\n")
    
    if(iter>1 && abs(loglik - loglik_trace[iter-1]) < tol) break
  }
  
  loglik_trace <- loglik_trace[1:iter]
  m <- G*(p+4) + (G-1)*q
  
  list(
    beta = beta,
    sigma2 = sigma2,
    lambda = lambda,
    nu = nu,
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

