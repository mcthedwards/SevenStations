# Function to sample signal parameters given noise parameters
gibbs_signal_alpha = function(yf, Xf, times.f, q, tau, beta) {
  
  # yf is original FTd data (i.e., signal + noise)
  # Xf is the FTd design matrix
  # D = 2 * pi * diag(f(lambda))
  # times.f is FTd time vector
  # beta is slope coefficient

  variances <- tau * q * 2 * pi  # 2 * pi * f, i.e., diagonal variance matrix
  invVar <- 1 / variances        # Inverse of diagonal variance matrix
  dummy <- matrix(NA, nrow = ncol(Xf), ncol = nrow(Xf)) 
  for (jj in 1:nrow(dummy)) {
    dummy[jj, ] <- invVar * Xf[, jj]
  }
  Sigma.inv <- dummy %*% Xf                          # t(X) %*% Dinv %*% X
  Sigma <- MASS::ginv(Sigma.inv)                     # Covariance matrix
  mu <- Sigma %*% dummy %*% (yf - beta * times.f)    # Mean vector
  alpha <- MASS::mvrnorm(1, mu, Sigma)               # Sample intercept(s) and slope
  
  return(alpha)
  
}

# Function to sample signal parameters given noise parameters
gibbs_signal_beta = function(yf, Xf, times.f, q, tau, alpha, 
                             beta.G, beta.v0) {
  
  # yf is original FTd data (i.e., signal + noise)
  # Xf is the FTd design matrix
  # D = 2 * pi * diag(f(lambda))
  # times.f is FTd time vector
  # alpha is y-intercepts for each changepoint
  # beta.G is global slope (i.e., mean of betas)
  # beta.v0 is variance of betas
  
  variances <- tau * q * 2 * pi  # 2 * pi * f, i.e., diagonal variance matrix
  invVar = 1 / variances
  inv.sigma2 = sum(times.f ^ 2 * invVar) + 1 / beta.v0
  sigma2 = 1 / inv.sigma2
  mu = (sum(times.f * (yf - Xf %*% alpha) * invVar) + beta.G) * sigma2
  
  # Dinv = diag(invVar)
  # a = t(times.f) %*% Dinv %*% times.f + 1 / beta.v0
  # noise.f = yf - Xf %*% alpha
  # mu = (t(times.f) %*% Dinv %*% noise.f + beta.G) / a
  # sigma2 = 1 / a
  
  beta = rnorm(1, mu, sqrt(sigma2))
  
  return(beta)
  
}




