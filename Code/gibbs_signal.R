# Function to sample signal parameters given noise parameters
gibbs_signal = function(yf, Xf, q, tau) {
  
  # yf is original FTd data (i.e., signal + noise)
  # Xf is the FTd design matrix
  # D = 2 * pi * diag(f(lambda))
  
  # To do: Cholesky or QR decomposition
  
  variances <- tau * q * 2 * pi  # 2 * pi * f, i.e., diagonal variance matrix
  invVar <- 1 / variances        # Inverse of diagonal variance matrix
  dummy <- matrix(NA, nrow = ncol(Xf), ncol = nrow(Xf)) 
  for (jj in 1:nrow(dummy)) {
    dummy[jj, ] <- invVar * Xf[, jj]
  }
  Sigma.inv <- dummy %*% Xf            # t(X) %*% Dinv %*% X
  Sigma <- MASS::ginv(Sigma.inv)       # Covariance matrix
  mu <- Sigma %*% dummy %*% yf         # Mean vector
  beta <- MASS::mvrnorm(1, mu, Sigma)  # Sample intercept(s) and slope
  
  return(beta)
  
  
}






