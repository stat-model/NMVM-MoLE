mix.reg.SMSN.EM <- function(y, X, R, G = 2,
                                 family = c("Skew.t","Skew.n","Skew.slash","Skew.cn"),
                                 max_iter = 200,
                                 tol = 1e-6,
                                 ridge = 1e-4,
                                 verbose = TRUE){
  
  family <- match.arg(family)
  n <- length(y)
  p <- ncol(X)
  q <- ncol(R)
  g <- G
  
  #------------------ Stable softmax (برداری) ------------------#
  softmax_mat <- function(M){
    M <- M - apply(M, 1, max)
    expM <- exp(M)
    expM / rowSums(expM)
  }
  
  #------------------ Log-densities (SMSN) ------------------#
  logSN <- function(e, sigma2, shape){
    z <- e / sqrt(sigma2)
    log(2) + dnorm(z, log=TRUE) + pnorm(shape*z, log.p=TRUE) - 0.5*log(sigma2)
  }
  
  logST <- function(e, sigma2, shape, nu){
    z <- e / sqrt(sigma2)
    log(2) + dt(z, df=nu, log=TRUE) +
      pt(shape*z*sqrt((nu+1)/(nu+z^2)), df=nu+1, log.p=TRUE) -
      0.5*log(sigma2)
  }
  
  logSlash <- function(e, sigma2, nu){
    z <- e^2 / sigma2
    log(nu) - (nu+1)*log(1+z/nu) - 0.5*log(sigma2)
  }
  
  logSCN <- function(e, sigma2, shape, nu){
    eps <- nu[1]; k <- nu[2]
    term1 <- log(eps) + dnorm(e, 0, sqrt(sigma2/k), log=TRUE) + pnorm(shape*e*sqrt(k/sigma2), log.p=TRUE)
    term2 <- log(1-eps) + dnorm(e, 0, sqrt(sigma2), log=TRUE) + pnorm(shape*e/sqrt(sigma2), log.p=TRUE)
    m <- pmax(term1, term2)
    m + log(exp(term1 - m) + exp(term2 - m))
  }
  
  logdens <- function(e, sigma2, shape, nu){
    switch(family,
           "Skew.n" = logSN(e,sigma2,shape),
           "Skew.t" = logST(e,sigma2,shape,nu),
           "Skew.slash" = logSlash(e,sigma2,nu),
           "Skew.cn" = logSCN(e,sigma2,shape,nu))
  }
  
  #------------------ Initialization ------------------#
  km <- kmeans(y, G)
  Class <- km$cluster
  
  beta <- matrix(0, p, g)
  sigma2 <- rep(var(y), g)
  shape <- rep(1, g)
  
  if(family=="Skew.t") nu <- rep(8,g)
  if(family=="Skew.slash") nu <- rep(5,g)
  if(family=="Skew.cn") nu <- matrix(c(0.9,10),2,g)
  if(family=="Skew.n") nu <- NULL
  
  alpha <- matrix(0, q, g-1)
  
  for(j in 1:g){
    idx <- which(Class==j)
    if(length(idx)>p){
      beta[,j] <- solve(crossprod(X[idx, ]) + ridge*diag(p),
                        crossprod(X[idx, ], y[idx]))
    }
  }
  
  loglik_trace <- numeric(max_iter)
  
  #================ EM =================#
  for(iter in 1:max_iter){
    
    #------------------ E-step برداری ------------------#
    mu_mat <- X %*% beta        # n x g
    dens_mat <- matrix(0, n, g)
    
    for(j in 1:g){
      if(family=="Skew.cn"){
        dens_mat[, j] <- exp(logSCN(y - mu_mat[,j], sigma2[j], shape[j], nu[,j]))
      } else if(family=="Skew.n"){
        dens_mat[, j] <- exp(logSN(y - mu_mat[,j], sigma2[j], shape[j]))
      } else if(family=="Skew.t"){
        dens_mat[, j] <- exp(logST(y - mu_mat[,j], sigma2[j], shape[j], nu[j]))
      } else if(family=="Skew.slash"){
        dens_mat[, j] <- exp(logSlash(y - mu_mat[,j], sigma2[j], nu[j]))
      }
    }
    
    pi_mat <- softmax_mat(cbind(R %*% alpha, 0))
    tau <- pi_mat * dens_mat
    tau <- pmax(tau, 1e-12)
    tau <- tau / rowSums(tau)
    
    #------------------ M-step ------------------#
    for(j in 1:g){
      w <- tau[, j]
      W <- diag(w)
      
      beta[,j] <- solve(t(X) %*% W %*% X + ridge*diag(p),
                        t(X) %*% W %*% y)
      
      e_j <- y - X %*% beta[,j]
      
      sigma2[j] <- sum(w * e_j^2) / sum(w)
      sigma2[j] <- max(sigma2[j], 1e-6)
      
      shape[j] <- sum(w * e_j^3)/sum(w*e_j^2)
      shape[j] <- max(min(shape[j], 10), -10)
      
      if(family %in% c("Skew.t","Skew.slash")){
        obj <- function(par) -sum(w * logdens(e_j, sigma2[j], shape[j], par))
        nu[j] <- optim(nu[j], obj, method="L-BFGS-B", lower=2.1, upper=50)$par
      }
      
      if(family=="Skew.cn"){
        obj <- function(par){
          eps <- plogis(par[1]); k <- exp(par[2])
          -sum(w * logSCN(e_j, sigma2[j], shape[j], c(eps,k)))
        }
        opt <- optim(c(qlogis(nu[1,j]), log(nu[2,j])), obj)
        nu[,j] <- c(plogis(opt$par[1]), exp(opt$par[2]))
      }
    }
    
    #------------------ M-step (gating) ------------------#
    loglik_alpha <- function(alpha_vec){
      alpha_mat <- matrix(alpha_vec, nrow=q, ncol=g-1)
      pi_mat <- softmax_mat(cbind(R %*% alpha_mat, 0))
      -sum(tau * log(pmax(pi_mat,1e-12)))
    }
    
    opt <- optim(as.vector(alpha), loglik_alpha, method="BFGS")
    alpha <- matrix(opt$par, nrow=q, ncol=g-1)
    
    #------------------ Log-likelihood برداری ------------------#
    loglik <- sum(log(rowSums(pi_mat * dens_mat)))
    loglik_trace[iter] <- loglik
    
    if(verbose) cat("Iter:", iter, "LogLik:", round(loglik,4), "\n")
    
    if(iter>1 && abs(loglik - loglik_trace[iter-1])<tol) break
  }
  
  loglik_trace <- loglik_trace[1:iter]
  
  #------------------ Model selection ------------------#
  d <- g*(p+2) + (g-1)*q
  if(family!="Skew.n") d <- d + g
  AIC <- -2*loglik + 2*d
  BIC <- -2*loglik + log(n)*d
  
  #------------------ Output ------------------#
  list(
    beta = beta,
    sigma2 = sigma2,
    shape = shape,
    nu = nu,
    alpha = alpha,
    tau = tau,
    loglik = loglik,
    loglik_trace = loglik_trace,
    AIC = AIC,
    BIC = BIC,
    EDC = -2*loglik + 0.2*sqrt(n)*d,
    CAIC = -2*loglik + 2*n*d/(n-d-1),
    ABIC = -2*loglik + d*log((n+2)/24),
    clusters = apply(tau,1,which.max),
    npar = d
  )
}

