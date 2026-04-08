fit_robust_bi_mole <- function(y, X, R, G=3, c=1.345, max_iter=200, tol=1e-6){
  
  if(!requireNamespace("nnet", quietly = TRUE)) stop("Package 'nnet' is required.")
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  
  # ---------- Initialization ----------
  beta <- matrix(0, p, G)
  sigma2 <- rep(var(y), G)
  alpha <- matrix(0, q, G-1)
  
  km <- kmeans(cbind(y, X), G)
  for(j in 1:G){
    idx <- which(km$cluster == j)
    if(length(idx) > p){
      XtX <- t(X[idx,]) %*% X[idx,]
      beta[,j] <- solve(XtX + 1e-6*diag(p)) %*% t(X[idx,]) %*% y[idx]
    }
  }
  
  softmax <- function(v){ exp(v - max(v)) / sum(exp(v - max(v))) }
  psi <- function(e){ ifelse(abs(e) <= c, e, c*sign(e)) }
  loglik_old <- -Inf
  
  # ---------- EM ----------
  for(iter in 1:max_iter){
    tau <- matrix(0, n, G)
    loglik <- 0
    
    # ===== E-step =====
    for(i in 1:n){
      eta_g <- c(R[i,] %*% alpha, 0)
      pi_i <- softmax(eta_g)
      tmp <- numeric(G)
      for(j in 1:G){
        e <- (y[i] - X[i,] %*% beta[,j]) / sqrt(sigma2[j])
        dens <- exp(-0.5 * psi(e)^2) / sqrt(sigma2[j])
        tmp[j] <- pi_i[j] * dens
      }
      den_i <- max(sum(tmp), 1e-12)
      tau[i,] <- tmp / den_i
      loglik <- loglik + log(den_i)
    }
    
    # ===== M-step: Experts =====
    for(j in 1:G){
      w <- tau[,j]
      for(k in 1:5){
        e <- y - X %*% beta[,j]
        w_rob <- psi(e) / pmax(abs(e), 1e-6)
        W <- w * w_rob
        XtWX <- t(X) %*% diag(c(W)) %*% X
        beta[,j] <- solve(XtWX + 1e-6*diag(p)) %*% t(X) %*% diag(c(W)) %*% y
      }
      e <- y - X %*% beta[,j]
      sigma2[j] <- max(sum(w*e^2)/sum(w), 1e-8)
    }
    
    # ===== M-step: Gating (multinom from nnet) =====
    if(G > 1){
      df <- data.frame(R)
      colnames(df) <- paste0("R", 1:ncol(R))
      df$cluster <- factor(apply(tau[,1:(G-1)],1,which.max))
      # Fit multinomial logit
      fit <- nnet::multinom(cluster ~ . , data=df, trace=FALSE)
      alpha <- t(coef(fit))  # q x (G-1)
    }
    
    # ===== Convergence =====
    if(abs(loglik - loglik_old) < tol) break
    loglik_old <- loglik
  }
  
  list(
    clusters = apply(tau,1,which.max),
    tau = tau,
    beta = beta,
    sigma2 = sigma2,
    alpha = alpha,
    loglik = loglik,
    npar = G*(p+1) + (G-1)*q,
    iter = iter
  )
}
