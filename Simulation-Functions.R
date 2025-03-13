# Simulation Function
rep.col <- function(x, n){
  matrix(rep(x, each = n), ncol = n, byrow = TRUE)
}

Coef.gen <- function(s, h, K, sig.beta, sig.eta, p, exact){
  beta0 <- c(rep(sig.beta, s), rep(0, p - s))
  theta <- rep.col(beta0, K)
  for(k in 1:K){
    if(exact){
      samp0 <- sample(1:p, h, replace = F)
      theta[samp0, k] <- theta[samp0, k] + rep(-sig.eta, h)
    }else{
      theta[1:100, k] <- theta[1:100, k] + rnorm(100, 0, h/100)
    }
  }
  return(list(theta = theta, beta0 = beta0))
}

ind.set <- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k == 1){
      ind.re <- c(ind.re, 1:n.vec[1])
    }else{
      ind.re <- c(ind.re, (sum(n.vec[1:(k - 1)]) + 1):sum(n.vec[1:k]))
    }
  }
  ind.re
}

library(mvtnorm)
Data.gen <- function(n0, nk, K, p, Coef.list, heter){
  B <- cbind(Coef.list$beta0, Coef.list$theta)
  n.vec <- c(n0, rep(nk, K))
  X <- NULL; y <- NULL; status <- NULL
  X <- rbind(X, rmvnorm(n.vec[1], rep(0, p), diag(1, p)))
  ind.k <- ind.set(n.vec, 1)
  logT <- X[ind.k, ] %*% B[, 1] + rnorm(n.vec[1])
  C <- runif(n.vec[1], 0, 12)
  delta <- (exp(logT) <= C) + 0
  y = c(y, logT*delta + log(C)*(1 - delta))
  status <- c(status, delta)
  if(heter == F){
    for(k in 2:(K + 1)){
      X <- rbind(X, rmvnorm(n.vec[k], rep(0, p), diag(1, p)))
      ind.k <- ind.set(n.vec, k)
      logT <- X[ind.k, ] %*% B[, k] + rnorm(n.vec[k])
      C <- runif(n.vec[k], 0, 12)
      delta <- (exp(logT) <= C) + 0
      y = c(y, logT*delta + log(C)*(1 - delta))
      status <- c(status, delta)
    }
  }else{
    for(k in 2:(K + 1)){
      Cov <- toeplitz(c(1, rep(1/k, 2*k - 3), rep(0, p - 2*(k - 1))))
      X <- rbind(X, rmvnorm(n.vec[k], rep(0, p), Cov))
      ind.k <- ind.set(n.vec, k)
      logT <- X[ind.k, ] %*% B[, k] + rnorm(n.vec[k])
      C <- runif(n.vec[k], 0, 12)
      delta <- (exp(logT) <= C) + 0
      y = c(y, logT*delta + log(C)*(1 - delta))
      status <- c(status, delta)
    }
  }
  return(list(X = X, y = y, status = status))
}

KM.Weight.Center <- function(n0, nk, K, Xylist){ 
  n.vec <- c(n0, rep(nk, K))
  X <- Xylist$X; y <- Xylist$y; status <- Xylist$status
  X.center <- NULL; y.center <- NULL
  for(k in 1:(K + 1)){
    ind.k <- ind.set(n.vec, k)
    y.k = y[ind.k]; X.k = X[ind.k, ]; status.k = status[ind.k]
    y.order = sort(y.k); X.order = X.k[order(y.k), ]; status.order = status.k[order(y.k)]
    w = NULL; n = length(y.order)
    w[1] = status.order[1]/n
    m = 1
    for(i in 2:n){
      m = m*((n - i + 1)/(n - i + 2))^status.order[i - 1]
      w[i] = status.order[i]/(n - i + 1)*m
    }
    y.bar <- sum(w*y.order)/sum(w); X.bar <- apply(w*X.order, 2, sum)/sum(w)
    X_minus <- X.order
    for(c in 1:n.vec[k]){
      X_minus[c, ] <- as.numeric(X_minus[c, ]) - as.numeric(X.bar)
    }
    X.center <- rbind(X.center, sqrt(n.vec[k]*w)*(X_minus))
    y.center <- c(y.center, sqrt(n.vec[k]*w)*(y.order - y.bar))
  }
  return(list(X = as.matrix(X.center), y = y.center))
}

library(glmnet)
Trans.AFT <- function(X, y, n0, nk, K){
  n.vec <- c(n0, rep(nk, K))
  p <- ncol(X)
  ind.1 <- 1:n.vec[1]
  cv.init <- cv.glmnet(X, y, nfolds = 8, lambda = seq(1, 0.1, length.out = 10)*sqrt(2*log(p)/sum(n.vec)))
  lam.const <- cv.init$lambda.min/sqrt(2*log(p)/sum(n.vec))
  theta.S <- as.numeric(glmnet(X, y, lambda = lam.const*sqrt(2*log(p)/sum(n.vec)))$beta) 
  theta.S <- theta.S*(abs(theta.S) >= lam.const*sqrt(2*log(p)/sum(n.vec)))
  # cv.eta <- cv.glmnet(X[ind.1, ], y[ind.1] - X[ind.1, ]%*%theta.S, lambda = seq(1, 0.1, length.out = 10)*sqrt(2*log(p)/length(ind.1)))
  # lam.const2 <- cv.eta$lambda.min/sqrt(2*log(p)/length(ind.1))
  # eta <- as.numeric(glmnet(x = X[ind.1,],y = y[ind.1] - X[ind.1,]%*%theta.S, lambda = lam.const2*sqrt(2*log(p)/length(ind.1)))$beta)
  eta <- as.numeric(glmnet(x = X[ind.1,],y = y[ind.1] - X[ind.1,]%*%theta.S, lambda = lam.const*sqrt(2*log(p)/length(ind.1)))$beta)
  eta <- eta*(abs(eta) >= lam.const*sqrt(2*log(p)/length(ind.1)))
  beta <- theta.S + eta
  list(beta = as.numeric(beta),theta.S = theta.S)
}

err.fun <- function(beta, est){
  est.err <- sum(abs(beta - est))
  #est.err <- sum((beta - est)^2)
  return(est.err)
}
