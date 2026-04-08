#========================================================
#   Unified Simulation Framework for NMVM & SMSN Models
#========================================================

#-----------------------------#
#   Utility: Stable Softmax   #
#-----------------------------#
softmax_stable <- function(eta){
  eta <- eta - max(eta)
  exp_eta <- exp(eta)
  exp_eta / sum(exp_eta)
}

#========================================================
#   Mixing Distributions (W)
#========================================================
mix.sample <- function(n, param, family = c("BS","Lindley","GIG","Exp")){
  
  family <- match.arg(family)
  
  if(family == "BS"){
    alpha <- param[1]
    U <- rnorm(n)
    W <- (alpha*U + sqrt((alpha*U)^2 + 4))^2 / 4
  }
  
  if(family == "Lindley"){
    alpha <- param[1]
    W <- VGAM::rlind(n, alpha)
  }
  
  if(family == "GIG"){
    W <- GIGrvg::rgig(n,
                      lambda = param[1],
                      chi    = param[2],
                      psi    = param[3])
  }
  
  if(family == "Exp"){
    W <- rexp(n, rate = 0.5)
  }
  
  return(W)
}

#========================================================
#   NMVM Generator
#========================================================
r.NMVM <- function(n, sigma2, lambda, theta,
                   family = c("Normal","SL","GHST","VG",
                              "NIG","GH","NMVBS","NMVL"),
                   skew = FALSE){
  
  family <- match.arg(family)
  
  # Base Gaussian
  Z <- rnorm(n, 0, sqrt(sigma2))
  
  # Mixing variable
  W <- switch(family,
              
              "Normal" = rep(1,n),
              
              "SL"     = mix.sample(n, NULL, "Exp"),
              
              "GHST"   = mix.sample(n, c(-theta/2, theta, 0), "GIG"),
              
              "VG"     = mix.sample(n, c(theta[1],0,theta[2]), "GIG"),
              
              "NIG"    = mix.sample(n, c(-0.5,1,theta^2), "GIG"),
              
              "GH"     = mix.sample(n, c(theta,theta,theta), "GIG"),
              
              "NMVBS"  = mix.sample(n, theta, "BS"),
              
              "NMVL"   = mix.sample(n, theta, "Lindley")
  )
  
  if(family == "Normal") lambda <- 0
  
  # Optional skew-normal innovation
  if(skew){
    Z <- sn::rsn(n, xi=0, omega=1, alpha=1)
  }
  
  Y <- lambda*W + sqrt(W)*Z
  
  return(Y)
}

#========================================================
#   NMVM Regression
#========================================================
r.NMVM.reg <- function(x, beta, sigma2, lambda, theta, family, skew=FALSE){
  eps <- r.NMVM(1, sigma2, lambda, theta, family, skew)
  as.numeric(x %*% beta + eps)
}

#========================================================
#   SMSN Generator (4 distributions)
#========================================================
r.SMSN.reg <- function(x, beta, sigma2, shape, nu,
                       family = c("ESN","ST","SSL","SCN")){
  
  family <- match.arg(family)
  mu <- as.numeric(x %*% beta)
  
  if(family == "ESN"){
    y <- sn::rsn(1, xi=mu, omega=sqrt(sigma2), alpha=shape, tau=nu)
  }
  
  if(family == "ST"){
    y <- sn::rst(1, xi=mu, omega=sqrt(sigma2), alpha=shape, nu=nu)
  }
  
  if(family == "SCN"){
    U <- ifelse(runif(1) < nu[1], nu[2], 1)
    z <- sn::rsn(1, xi=0, omega=1, alpha=shape)
    y <- mu + sqrt(sigma2/U)*z
  }
  
  if(family == "SSL"){
    U <- rbeta(1, nu, 1)
    z <- sn::rsn(1, xi=0, omega=1, alpha=shape)
    y <- mu + sqrt(sigma2/U)*z
  }
  
  return(y)
}

#========================================================
#   Mixture of Experts: NMVM
#========================================================
r.mix.NMVM <- function(n, X, beta, sigma2, lambda, theta,
                       alpha, R, family, skew=FALSE){
  
  g <- length(sigma2)
  y <- numeric(n)
  z <- numeric(n)
  
  for(i in 1:n){
    
    eta <- c(R[i,] %*% alpha, 0)
    pi_i <- softmax_stable(eta)
    
    comp <- sample(1:g, size=1, prob=pi_i)
    z[i] <- comp
    
    y[i] <- r.NMVM.reg(X[i,], beta[,comp],
                       sigma2[comp], lambda[comp],
                       theta[,comp], family, skew)
  }
  
  list(y=y, class=z)
}

#========================================================
#   Mixture of Experts: SMSN
#========================================================
r.mix.SMSN <- function(n, X, beta, sigma2, shape, nu,
                       alpha, R, family){
  
  g <- length(sigma2)
  y <- numeric(n)
  z <- numeric(n)
  
  for(i in 1:n){
    
    eta <- c(R[i,] %*% alpha, 0)
    pi_i <- softmax_stable(eta)
    
    comp <- sample(1:g, 1, prob=pi_i)
    z[i] <- comp
    
    y[i] <- r.SMSN.reg(X[i,], beta[,comp],
                       sigma2[comp], shape[comp],
                       nu[,comp], family)
  }
  
  list(y=y, class=z)
}

#========================================================
#   Prediction (Expectation-based)
#========================================================
predict.mix.NMVM <- function(X, beta, sigma2, lambda, theta,
                             alpha, R, family){
  
  n <- nrow(X)
  g <- length(sigma2)
  
  yhat <- numeric(n)
  
  for(i in 1:n){
    eta <- c(R[i,] %*% alpha, 0)
    pi_i <- softmax_stable(eta)
    
    mu_j <- sapply(1:g, function(j)
      X[i,] %*% beta[,j] + lambda[j])
    
    yhat[i] <- sum(pi_i * mu_j)
  }
  
  return(yhat)
}

#========================================================
#   Matrix Inversion (Stable)
#========================================================
inv.mat <- function(M){
  eig <- eigen(M)
  eig$vectors %*% diag(1/eig$values) %*% t(eig$vectors)
}