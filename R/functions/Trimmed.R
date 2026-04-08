fit_trimmed_mole <- function(y, X, R, g=3, alpha_trim=0.1,
                            max_iter=200, tol=1e-6){
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  h <- floor((1 - alpha_trim)*n)
  
  # ---------- Initialization ----------
  km <- kmeans(cbind(y,X), g)
  z <- km$cluster
  
  beta <- matrix(0, p, g)
  sigma2 <- rep(var(y), g)
  alpha_g <- matrix(0, q, g-1)
  
  for (j in 1:g) {
    idx <- which(z == j)
    if (length(idx) > p) {
      XtX <- t(X[idx, ]) %*% X[idx, ]
      beta[, j] <- solve(XtX + 1e-6*diag(p)) %*% t(X[idx, ]) %*% y[idx]
    }
  }
  
  softmax <- function(v){
    exp(v - max(v)) / sum(exp(v - max(v)))
  }
  
  loglik_old <- -Inf
  
  # ---------- EM ----------
  for(iter in 1:max_iter){
    
    log_comp <- matrix(0, n, g)
    
    # ===== E-step =====
    for(i in 1:n){
      
      eta <- c(R[i,] %*% alpha_g, 0)
      pi_i <- softmax(eta)
      
      for(j in 1:g){
        mu <- X[i,] %*% beta[,j]
        log_comp[i,j] <- log(pi_i[j] + 1e-12) +
          dnorm(y[i], mu, sqrt(sigma2[j]), log=TRUE)
      }
    }
    
    z_new <- apply(log_comp,1,which.max)
    score <- apply(log_comp,1,max)
    
    # ===== Trimming =====
    keep <- order(score, decreasing=TRUE)[1:h]
    z_trim <- z_new
    z_trim[-keep] <- 0
    
    # ===== M-step: Experts =====
    for(j in 1:g){
      idx <- which(z_trim == j)
      
      if(length(idx) > p){
        Xj <- X[idx,]
        yj <- y[idx]
        
        XtX <- t(Xj) %*% Xj
        beta[,j] <- solve(XtX + 1e-6*diag(p)) %*% t(Xj) %*% yj
        
        e <- yj - Xj %*% beta[,j]
        sigma2[j] <- max(mean(e^2), 1e-8)
      }
    }
    
    # ===== M-step: Gating =====
    # فقط روی داده‌های نگه‌داشته‌شده
    for(j in 1:(g-1)){
      idx <- keep
      z_logit <- as.numeric(z_trim[idx] == j)
      
      fit <- glm(z_logit ~ R[idx,] - 1, family=binomial())
      alpha_g[,j] <- coef(fit)
    }
    
    # ===== Log-likelihood =====
    loglik <- sum(score[keep])
    
    if(abs(loglik - loglik_old) < tol) break
    loglik_old <- loglik
  }
  
  list(
    clusters = z_trim,
    beta = beta,
    sigma2 = sigma2,
    alpha = alpha_g,
    loglik = loglik,
    npar = g*(p+1) + (g-1)*q,
    iter = iter
  )
}