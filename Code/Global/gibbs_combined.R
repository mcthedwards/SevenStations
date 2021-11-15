

gibbs_temperature_overall = function(data,
                             Ntotal,
                             burnin,
                             thin = 1,
                             M = 1,
                             g0.alpha = 1,
                             g0.beta = 1,
                             k.theta = 0.01,
                             tau.alpha = 0.001,
                             tau.beta = 0.001,
                             mu.0 = 0,
                             sigma2.0 = 0.1,
                             a0 = 0.1,  # Vague enough?
                             b0 = 0.1,  # Vague enough?
                             kmax = 100,
                             L = 20,        #L = max(20, length(data) ^ (1 / 3)),   have just changed to 20 as lengths differ. <20 anyway
                             printerval = 1000
                             #,
                             #mu_0 = 0.01,
                             #sigma_0 = 0.001
                            ) {
  
  
  #####
  # Preamble
  #####
  
  
  if (burnin >= Ntotal) stop("burnin must be less than Ntotal")
  if (any(c(M, g0.alpha, g0.beta, tau.alpha, tau.beta, k.theta) <= 0)) stop("M, g0.alpha, g0.beta, tau.alpha, tau.beta, and k.theta must be strictly positive")
  if (any(c(Ntotal, thin, kmax, L) %% 1 != 0) || any(c(Ntotal, thin, kmax, L) <= 0)) stop("Ntotal, thin, kmax, and L must be strictly positive integers")
  if ((burnin %% 1 != 0) || (burnin < 0)) stop("burnin must be a non-negative integer")
  
  
  # Open objects for storage
  # need noise parameters for each location. Each row is a diff location.
  tau <- matrix(NA, nrow = 7, ncol = Ntotal)
  k <- matrix(NA, nrow = 7, ncol = Ntotal)
  
  # store v and W as list of matrices
  V <- lapply(1:7, function(x) matrix(NA, nrow = L, ncol = Ntotal))
  W <- lapply(1:7, function(x) matrix(NA, nrow = L + 1, ncol = Ntotal))
  
  # store beta in list (one beta for each location). Each component of the list will be a vector.
  betas <- vector(mode = "list", length = 7) 
  
  #intercepts for each location. each component of the list will be a matrix.
  alphas <- vector(mode = "list", length = 7)
  
  #global slope vector
  beta.G <- rep(NA, Ntotal) 
  beta.v0 <- rep(NA, Ntotal)
  
  #set initial value of beta.G
  beta.G[1] = 0.01
  
  beta.v0[1] = 1 / rgamma(1, a0, b0)  # Initial value from hyperprior
  
  # matrices to store Xt_i, Xf_i and yf_i
  Xts <- vector(mode = "list", length = 7)
  Xfs <- vector(mode = "list", length = 7)
  yfs <- vector(mode = "list", length = 7)
  
  # lists for proposal parameters
  eps_list <- vector(mode = "list", length = 7) 
  omega_list <- vector(mode = "list", length = 7)  
  lambda_list <- vector(mode = "list", length = 7) 
  
  
  #beta densities. This will become a list of lists as beta densities for each location will be a list itself
  dbs_list <- vector(mode = "list", length = 7) 
  
  #split data into separate dfs for each region
  data_split = split(data, data$location, drop = T)
  
  #vector for location means which we use at the end to compute the psds 
  location_means <- rep(NA, 7)
  
  #times.f for each location
  times.f_list <- vector(mode = "list", length = 7) 
  
  
  #location loop to pre-process, set intitial values etc.
  
  for (i in 1:7) {
    
    n = nrow(data_split[[i]])
    
    # Which boundary frequencies to remove from likelihood computation and tau sample
    if (n %% 2) {  # Odd length time series
      bFreq = 1  # Remove first 
    } 
    else {  # Even length time series
      bFreq = c(1, n)  # Remove first and last
    }
    
    
    times = data_split[[i]]$year
    temp = data_split[[i]]$temp
    
    # Change points
    cp = c()
    cp[1] = 1
    for (j in 2:nrow(data_split[[i]])) {
      cp[j] = (data_split[[i]]$site_number[j] != data_split[[i]]$site_number[j - 1])
    }
    cp = which(cp == 1)
    
    
    # Define changepoints as where the new level starts
    # changepoints = c(1, changepoints) 
    cpn = length(cp)
    
    # Time domain design matrix
    Xt = matrix(0, nrow = length(times), ncol = cpn)
    for (j in 1:(cpn-1)) {
      Xt[cp[j]:(cp[j+1]-1), j] = 1  # Indicators created to show which part of intercept to estimate
    }
    Xt[cp[cpn]:length(times), cpn] = 1  # Populating last segment with 1s
    
    # Note: Keep slope detached from design matrix
    
    # Frequency domain design matrix
    Xf = matrix(NA, nrow = length(times), ncol = cpn )
    for (j in 1:(cpn)) {
      Xf[, j] = fast_ft(Xt[, j])
    }
    
    # Fourier transform time explanatory variable
    times.mean = mean(times)
    times = times - mean(times)  # Mean center times (numerical stability)
    data_split[[i]]$times = times
    
    times.f = fast_ft(times)  
    
    times.f_list[[i]] = times.f
    
    # Original data in frequency domain
    temp.mean = mean(temp)
    location_means[i] = temp.mean
    temp = temp - mean(temp)  # Mean center
    data_split[[i]]$temp = temp
    
    # Original data in frequency domain
    yf = fast_ft(temp) 
    
    
    # Store Xf, Xt and yf for each location
    Xts[[i]] = Xt
    Xfs[[i]] = Xf
    yfs[[i]] = yf
    
    omega <- 2 * (1:(n / 2 + 1) - 1) / n  # Frequencies on unit interval
    omega_list[[i]] <- omega
    
    
    # do we need lambda. not used anywhere. just omega.
    lambda <- pi * omega  # Angular frequencies on [0, pi]
    

    # Noise parameters
    # Starting values
    tau[i,1] <- stats::var(temp) / (2 * pi)
    k[i,1] = sample(1:kmax, 1)  
    
    # Optimise starting values for DP parameters
    V[[i]][, 1] = vFromP(rep(1 / (L + 1), L))
    W[[i]][, 1] = seq(from= 1 / (2 * k[1]), to = 1 - 1 / (2 * k[1]), length.out = L + 1)
    
    
    #beta
    betas[[i]] = rep(NA, Ntotal)
    
    #alpha (intercepts) matrix for each location - depends on the number of change points
    alphas[[i]] <- matrix(NA, nrow = Ntotal, ncol = cpn)

    
    # Use least squares estimate as starting point for signal
    ##lm.fit = lm(temp ~ Xt - 1)
    ##betas[[i]][1, ] = as.numeric(coef(lm.fit))
    
    # Use least squares estimate as starting point for signal
    lm.fit = lm(temp ~ cbind(Xt, times) - 1)
    alphas[[i]][1, ] = as.numeric(coef(lm.fit))[1:cpn]
    betas[[i]][1] = as.numeric(coef(lm.fit))[cpn + 1]
    
    
    # Metropolis proposal parameters for V, U, W, Z.
    eps <- seq(1, L + 1) / (seq(1, L + 1) + 2 * sqrt(n))  
    eps_list[[i]] <- eps
    
    
    
    #beta densities - list of lists 
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
    
    
    # each location has different beta densities
    dbs_list[[i]] = db.list
    
    
    }
  }
  
    
  #loop over steps
  #there will be a location loop nested inside this for signal and noise, and then a global parameters block outside of the locations loop.
  for (i in 2:Ntotal) {
      
    if (i %% printerval == 0) print(paste0("Iteration ", i, " of ", Ntotal))
      
      
    for (j in 1:7) {
      
        
      #####
      # 1. Noise block
      #####
      
      # Construct noise series by subtracting signal series
      #noise <- data_split[[j]]$temp - as.vector(Xts[[j]] %*% betas[[j]][i - 1, ]) 
      signal = as.vector(Xts[[j]] %*% alphas[[j]][i - 1, ]) + betas[[j]][i - 1] * data_split[[j]]$times
      noise <- data_split[[j]]$temp - signal
      FZ <- fast_ft(noise)  # Frequency domain
      pdgrm <- abs(FZ) ^ 2  # Periodogram
      
      # Sample noise parameters
      # add j indices to this bit 
      theta = gibbs_noise(omega_list[[j]], FZ, k[j, i - 1], V[[j]][,i - 1], W[[j]][,i - 1], tau[j,i - 1], L, 
                          M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                          kmax, pdgrm, dbs_list[[j]], eps_list[[j]])
      
      # Store noise parameters
      k[j, i] = theta$k
      V[[j]][, i] = theta$V
      W[[j]][, i] = theta$W
      tau[j, i] = theta$tau
      q = theta$q                   #q is just an input into the signal block? 
      
    
      
      #####
      # 2. Signal block
      #####
      
      # Sample signal parameters
      
      # (a). y-intercepts (assumes U(-Inf, Inf) prior)
      alphas[[j]][i, ] <- gibbs_signal_alpha(yfs[[j]], Xfs[[j]], times.f_list[[j]], q, tau[j,i], betas[[j]][i - 1])
      
      # (b). Slope (assumes N(beta.G, beta.v0) prior)
      betas[[j]][i] <- gibbs_signal_beta(yfs[[j]], Xfs[[j]], times.f_list[[j]], q, tau[j,i], alphas[[j]][i, ], 
                                   beta.G[i - 1], beta.v0[i - 1])
      
      
      
  }
    #####
    # 3. Global parameters block
    #####
    beta.i = sapply(betas, function(x) x[i])
    beta.G[i] = gibbs_global_slope(beta.i, beta.v0[i-1], mu.0, sigma2.0)
    
    beta.v0[i] = gibbs_global_variance(beta.i, beta.G[i], a0, b0)
    
  }
    
    
  # Which iterations to keep
  keep <- seq(burnin + 1, Ntotal, by = thin)
  k <- k[,keep]
  tau <- tau[,keep]
  V <- lapply(V, function(x) x[,keep])
  W <- lapply(W, function(x) x[,keep])
  betas <- lapply(betas, function(x) x[keep])
  beta.G <- beta.G[keep]
  beta.v0 <- beta.v0[keep]
  
  
  #now loop over seven locations again for psds, recons
  psd_list <- vector(mode = "list", length = 7)
  recon_list <- vector(mode = "list", length = 7)
  
  for (i in 1:7){
  
    #each location has different omega
    fpsd.sample <- matrix(NA, nrow = length(omega_list[[i]]), ncol = length(keep))
    recon <- matrix(NA, nrow = nrow(data_split[[i]]), ncol = length(keep))
    
    # Store PSDs and reconstructed signals
    for (isample in 1:length(keep)) {
      fpsd.sample[, isample] <- tau[isample] * qpsd(omega_list[[i]], 
                                                    V[[i]][, isample], 
                                                    W[[i]][, isample], 
                                                    k[isample], 
                                                    dbs_list[[i]])$psd
      recon[, isample] <- as.vector(Xts[[i]] %*% alphas[[i]][isample, ]) +
        betas[[i]][isample] * data_split[[i]]$times + 
        location_means[i]
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
    
    psd_list[[i]] = list(psd.median = psd.median,
                         psd.mean = psd.mean,
                         psd.p05 = psd.p05,
                         psd.p95 = psd.p95)
    
    recon_list[[i]] = list(recon.median = recon.median,
                           recon.mean = recon.mean,
                           recon.p05 = recon.p05,
                           recon.p95 = recon.p95)
  }
  
  # List to output
  output = list(psd_list = psd_list,
                recon_list = recon_list,
                k = k,
                tau = tau,
                V = V,
                W = W, 
                alphas = alphas,
                betas = betas,
                beta.G = beta.G,
                beta.v0 = beta.v0)
  
  class(output) = "temp"  # Assign S3 class to object
  
  return(output)  # Return output
    
  }
  
  

  