
gibbs_psd = function(data,
                     Ntotal,
                     burnin,
                     thin = 1,
                     M = 1,
                     g0.alpha = 1,
                     g0.beta = 1,
                     k.theta = 0.01,
                     tau.alpha = 0.001,
                     tau.beta = 0.001,
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
  
  # Mean center
  if (abs(mean(data)) > tol) {
    data <- data - mean(data)
    # warning("data has been mean-centered")
  }
  
  # Optimal rescaling to prevent numerical issues
  rescale = stats::sd(data)
  data = data / rescale  # Data now has standard deviation 1
  
  if (burnin >= Ntotal) stop("burnin must be less than Ntotal")
  if (any(c(M, g0.alpha, g0.beta, tau.alpha, tau.beta, k.theta) <= 0)) stop("M, g0.alpha, g0.beta, tau.alpha, tau.beta, and k.theta must be strictly positive")
  if (any(c(Ntotal, thin, kmax, L) %% 1 != 0) || any(c(Ntotal, thin, kmax, L) <= 0)) stop("Ntotal, thin, kmax, and L must be strictly positive integers")
  if ((burnin %% 1 != 0) || (burnin < 0)) stop("burnin must be a non-negative integer")
  
  FZ = fast_ft(data)  # FFT data to frequency domain.  NOTE: Must be mean-centered.
  
  pdgrm = abs(FZ) ^ 2   # Periodogram: NOTE: the length is n here.
  
  omega <- 2 * (1:(n / 2 + 1) - 1) / n  # Frequencies on unit interval
  lambda <- pi * omega  # Angular frequencies on [0, pi]
  
  # Open objects for storage
  tau <- rep(NA, Ntotal)
  V <- matrix(NA, nrow = L, ncol = Ntotal)
  W <- matrix(NA, nrow = L + 1, ncol = Ntotal)
  k <- rep(NA, Ntotal)
  
  # Starting values
  tau[1] <- stats::var(data) / (2 * pi)
  k[1] = sample(1:kmax, 1)  
  
  # Optimise starting values for DP parameters
  V[, 1] = vFromP(rep(1 / (L + 1), L))
  W[, 1] = seq(from=1 / (2 * k[1]), to = 1 - 1 / (2 * k[1]), length.out = L + 1)
  
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
    
    # Sample noise parameters
    theta = gibbs_noise(omega, FZ, k[i - 1], V[, i - 1], W[, i - 1], tau[i - 1], L, 
                        M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                        kmax, pdgrm, db.list, eps)
    
    # Store noise parameters
    k[i] = theta$k
    V[, i] = theta$V
    W[, i] = theta$W
    tau[i] = theta$tau
    
  }
  
  # Which iterations to keep
  keep <- seq(burnin + 1, Ntotal, by = thin)
  k <- k[keep]
  tau <- tau[keep]
  V <- V[, keep]
  W <- W[, keep]
  
  fpsd.sample <- matrix(NA, nrow = length(omega), ncol = length(keep))
  
  # Store PSDs
  for (isample in 1:length(keep)) {
    fpsd.sample[, isample] <- tau[isample] * qpsd(omega, 
                                                  V[, isample], 
                                                  W[, isample], 
                                                  k[isample], 
                                                  db.list)$psd
  }
  
  # Compute point estimates and 90% Pointwise CIs
  psd.median <- apply(fpsd.sample, 1, stats::median)
  psd.mean <- apply(fpsd.sample, 1, mean)
  psd.p05 <- apply(fpsd.sample, 1, stats::quantile, probs=0.05)
  psd.p95 <- apply(fpsd.sample, 1, stats::quantile, probs=0.95)
  
  # Compute periodogram
  N = length(psd.median)  # N = (n + 1) / 2 (ODD) or N = n / 2 + 1 (EVEN)
  pdgrm = (abs(stats::fft(data)) ^ 2 / (2 * pi * n))[1:N]
  
  # List to output
  output = list(psd.median = psd.median * rescale ^ 2,
                psd.mean = psd.mean * rescale ^ 2,
                psd.p05 = psd.p05 * rescale ^ 2,
                psd.p95 = psd.p95 * rescale ^ 2,
                k = k,
                tau = tau,
                V = V,
                W = W, 
                pdgrm = pdgrm * rescale ^ 2,
                n = n)
  
  class(output) = "psd"  # Assign S3 class to object
  
  return(output)  # Return output
  
}
