fit_contaminated_mole <- function(y, X, R, G=3,
                                           epsilon=rep(0.05,G),
                                           eta=rep(10,G),
                                           max_iter=200, tol=1e-6){
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  
  beta <- matrix(0,p,G)
  sigma2 <- rep(var(y), G)
  alpha <- matrix(0,q,G-1)
  
  # Init
  km <- kmeans(cbind(y,X), G)
  for (j in 1:G){
    idx <- which(km$cluster==j)
    if(length(idx) > p){
      XtX <- t(X[idx,]) %*% X[idx,]
      beta[,j] <- solve(XtX + 1e-6*diag(p)) %*% t(X[idx,]) %*% y[idx]
    }
  }
  
  softmax <- function(v){
    exp(v - max(v)) / sum(exp(v - max(v)))
  }
  
  loglik_old <- -Inf
  
  for(iter in 1:max_iter){
    
    tau <- matrix(0,n,G)
    v <- matrix(0,n,G)
    loglik <- 0
    
    # ===== E-step =====
    for(i in 1:n){
      eta_g <- c(R[i,] %*% alpha, 0)
      pi_i <- softmax(eta_g)
      
      tmp <- numeric(G)
      for(j in 1:G){
        mu <- X[i,] %*% beta[,j]
        e <- y[i] - mu
        f1 <- dnorm(e, 0, sqrt(sigma2[j]))
        f2 <- dnorm(e, 0, sqrt(eta[j]*sigma2[j]))
        denom <- max((1-epsilon[j])*f1 + epsilon[j]*f2, 1e-12)
        v[i,j] <- epsilon[j]*f2 / denom
        tmp[j] <- pi_i[j] * denom
      }
      tau[i,] <- tmp / sum(tmp)
      loglik <- loglik + log(sum(tmp))
    }
    
    # ===== M-step: Experts =====
    for(j in 1:G){
      w <- tau[,j]*(1 - v[,j] + v[,j]/eta[j])
      beta[,j] <- solve(t(X) %*% (w*X) + 1e-6*diag(p)) %*% t(X) %*% (w*y)
      e <- y - X %*% beta[,j]
      sigma2[j] <- max(sum(w*e^2)/sum(tau[,j]), 1e-8)
    }
    
    # ===== M-step: Gating =====
    for(j in 1:(G-1)){
      z <- tau[,j]
      fit <- glm(z ~ R - 1, family=binomial())
      alpha[,j] <- coef(fit)
    }
    
    # ===== Convergence =====
    if(abs(loglik-loglik_old) < tol) break
    loglik_old <- loglik
  }
  
  list(
    clusters = apply(tau,1,which.max),
    tau = tau,
    beta = beta,
    sigma2 = sigma2,
    alpha = alpha,
    loglik = loglik,
    npar = G*(p+2) + (G-1)*q,
    iter = iter
  )
}
