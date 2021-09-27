
gibbs_temperature = function(data,
                             times,  
                             changepoints,
                             Ntotal,
                             burnin,
                             thin = 1,
                             M = 1,
                             g0.alpha = 1,
                             g0.beta = 1,
                             k.theta = 0.01,
                             tau.alpha = 0.001,
                             tau.beta = 0.001,
                             beta.v0 = 1,  
                             mu.0 = 0,
                             sigma2.0 = 1,
                             kmax = 100,
                             L = max(20, length(data) ^ (1 / 3)),
                             printerval = 1000) {
  
  #####
  # Preamble
  #####
  
  n = length(data)
  
  # Which boundary frequencies to remove from likelihood computation and tau sample
  if (n %% 2) {  # Odd length time series
    bFreq = 1  # Remove first 
  } 
  else {  # Even length time series
    bFreq = c(1, n)  # Remove first and last
  }
  
  # Tolerance for mean centering
  tol <- 1e-4
  
  if (burnin >= Ntotal) stop("burnin must be less than Ntotal")
  if (any(c(M, g0.alpha, g0.beta, tau.alpha, tau.beta, k.theta) <= 0)) stop("M, g0.alpha, g0.beta, tau.alpha, tau.beta, and k.theta must be strictly positive")
  if (any(c(Ntotal, thin, kmax, L) %% 1 != 0) || any(c(Ntotal, thin, kmax, L) <= 0)) stop("Ntotal, thin, kmax, and L must be strictly positive integers")
  if ((burnin %% 1 != 0) || (burnin < 0)) stop("burnin must be a non-negative integer")
  
  # Define changepoints as where the new level starts
  cpn = length(changepoints)
  
  # Time domain design matrix
  Xt = matrix(0, nrow = length(times), ncol = cpn)
  for (j in 1:(cpn-1)) {
    Xt[changepoints[j]:(changepoints[j+1]-1), j] = 1  # Indicators created to show which part of intercept to estimate
  }
  Xt[changepoints[cpn]:length(times), cpn] = 1  # Populating last segment with 1s
  # Note: Keep slope detached from design matrix

  # Frequency domain design matrix
  Xf = matrix(NA, nrow = length(times), ncol = cpn)
  for (j in 1:(cpn)) {
    Xf[, j] = fast_ft(Xt[, j])
  }
  
  # Fourier transform time explanatory variable
  times.mean = mean(times)
  times = times - mean(times)  # Mean center times (numerical stability)
  times.f = fast_ft(times)  
  
  # Original data in frequency domain
  data.mean = mean(data)
  data = data - mean(data)  # Mean center
  yf = fast_ft(data) 
  
  omega <- 2 * (1:(n / 2 + 1) - 1) / n  # Frequencies on unit interval
  lambda <- pi * omega  # Angular frequencies on [0, pi]
  
  # Open objects for storage
  tau <- rep(NA, Ntotal)
  V <- matrix(NA, nrow = L, ncol = Ntotal)
  W <- matrix(NA, nrow = L + 1, ncol = Ntotal)
  k <- rep(NA, Ntotal)
  
  alpha <- matrix(NA, nrow = Ntotal, ncol = cpn)
  beta <- rep(NA, Ntotal)
  beta.G <- rep(NA, Ntotal)
  
  # Starting values
  tau[1] <- stats::var(data) / (2 * pi)
  k[1] = sample(1:kmax, 1)  
  
  # Optimise starting values for DP parameters
  V[, 1] = vFromP(rep(1 / (L + 1), L))
  W[, 1] = seq(from=1 / (2 * k[1]), to = 1 - 1 / (2 * k[1]), length.out = L + 1)
  
  # Use least squares estimate as starting point for signal
  lm.fit = lm(data ~ cbind(Xt, times) - 1)
  alpha[1, ] = as.numeric(coef(lm.fit))[1:cpn]
  beta[1] = beta.G[1] = as.numeric(coef(lm.fit))[cpn + 1]
  
  # Metropolis proposal parameters for V, U, W, Z.
  eps <- seq(1, L + 1) / (seq(1, L + 1) + 2 * sqrt(n))  
  
  # Pre-populate beta densities up until kmax
  db.names <- paste("db", 1:kmax, sep = "")
  db.list <- vector("list", kmax)
  names(db.list) <- db.names
  for (kk in 1:kmax) {
    db.list[[kk]] <- matrix(dbeta(omega,
                                  rep(1:kk, each = length(omega)),
                                  rep(kk:1, each = length(omega))),
                            ncol = length(omega),
                            byrow = TRUE)
  }
  
  for (i in 2:Ntotal) {
    
    if (i %% printerval == 0) print(paste0("Iteration ", i, " of ", Ntotal))
    
    #####
    # 1. Noise block
    #####
    
    # Construct noise series by subtracting signal series
    signal = as.vector(Xt %*% alpha[i - 1, ]) + beta[i - 1] * times
    noise <- data - signal
    FZ <- fast_ft(noise)  # Frequency domain
    pdgrm <- abs(FZ) ^ 2  # Periodogram
    
    # Sample noise parameters
    theta = gibbs_noise(omega, FZ, k[i - 1], V[, i - 1], W[, i - 1], tau[i - 1], L, 
                        M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                        kmax, pdgrm, db.list, eps)
    
    # Store noise parameters
    k[i] = theta$k
    V[, i] = theta$V
    W[, i] = theta$W
    tau[i] = theta$tau
    q = theta$q
    
    #####
    # 2. Signal block
    #####
    
    # Sample signal parameters
    
    # (a). y-intercepts (assumes U(-Inf, Inf) prior)
    alpha[i, ] <- gibbs_signal_alpha(yf, Xf, times.f, q, tau[i], beta[i - 1])
    
    # (b). Slope (assumes N(beta.G, beta.v0) prior)
    beta[i] <- gibbs_signal_beta(yf, Xf, times.f, q, tau[i], alpha[i, ], 
                                 beta.G[i - 1], beta.v0)
    
    #####
    # 3. Global parameters block
    #####
    
    beta.G[i] = gibbs_global(beta[i], beta.v0, mu.0, sigma2.0)
    # To do: tweak this such that beta.G takes account of all stations
    # I.e., beta[i] should be a vector of length 7 rather than 1 here
    # And it should be outside the local loops
    
  }
  
  # Which iterations to keep
  keep <- seq(burnin + 1, Ntotal, by = thin)
  k <- k[keep]
  tau <- tau[keep]
  V <- V[, keep]
  W <- W[, keep]
  beta <- beta[keep]
  beta.G <- beta.G[keep]
  
  fpsd.sample <- matrix(NA, nrow = length(omega), ncol = length(keep))
  recon <- matrix(NA, nrow = n, ncol = length(keep))
  
  # Store PSDs and reconstructed signals
  for (isample in 1:length(keep)) {
    fpsd.sample[, isample] <- tau[isample] * qpsd(omega, 
                                                  V[, isample], 
                                                  W[, isample], 
                                                  k[isample], 
                                                  db.list)$psd
    recon[, isample] <- as.vector(Xt %*% alpha[isample, ]) +
      beta[isample] * times + 
      data.mean
  }
  
  # Compute point estimates and 90% Pointwise CIs
  psd.median <- apply(fpsd.sample, 1, stats::median)
  psd.mean <- apply(fpsd.sample, 1, mean)
  psd.p05 <- apply(fpsd.sample, 1, stats::quantile, probs=0.05)
  psd.p95 <- apply(fpsd.sample, 1, stats::quantile, probs=0.95)
  
  recon.median <- apply(recon, 1, stats::median)
  recon.mean <- apply(recon, 1, mean)
  recon.p05 <- apply(recon, 1, stats::quantile, probs=0.05)
  recon.p95 <- apply(recon, 1, stats::quantile, probs=0.95)

  # List to output
  output = list(psd.median = psd.median,
                psd.mean = psd.mean,
                psd.p05 = psd.p05,
                psd.p95 = psd.p95,
                recon.median = recon.median,
                recon.mean = recon.mean,
                recon.p05 = recon.p05,
                recon.p95 = recon.p95,
                k = k,
                tau = tau,
                V = V,
                W = W, 
                alpha = alpha,
                beta = beta,
                beta.G = beta.G,
                data = data + data.mean,
                times = times + times.mean)
  
  class(output) = "temp"  # Assign S3 class to object
  
  return(output)  # Return output
  
}
