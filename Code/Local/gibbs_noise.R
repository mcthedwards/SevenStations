# Function to sample noise parameters given signal parameters
gibbs_noise = function(omega, FZ, k, V, W, tau, L, 
                       M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                       kmax, pdgrm, db.list, eps) {
  
  n = length(FZ)
  
  f.store <- lpost(omega, FZ, V, W, k, tau,
                   M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                   pdgrm, db.list)   
  
  #####
  # Step 1: Metropolis proposal for k
  #####
  bold <- stats::runif(1)
  if (bold < 0.75) {  # 75/25 split
    jump <- sample(-1:1, 1, prob = rep(1 / 3, 3))  # Ordinary proposal
  }
  else {
    jump <- round(stats::rt(1, 1))  # Bold proposal - discrete Cauchy
  }
  k.star <- k + jump
  while (k.star < 1 || k.star > kmax) {  # A bit hacky to ensure k doesn't go out of bounds
    if (bold < 0.75) {
      jump <- sample(-1:1, 1, prob = rep(1 / 3, 3))  # Ordinary proposal
    }
    else {
      jump <- round(stats::rt(1, 1))  # Bold proposal
    }
    k.star <- k + jump
  }
  
  # log posterior for proposal
  f.k.star <- lpost(omega, FZ, V, W, k.star, tau,
                    M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                    pdgrm, db.list)   
  
  # log posterior of previous iteration
  f.k <- f.store
  
  #####
  # Accept/reject
  alpha1 <- min(0, f.k.star - f.k)  # log acceptance ratio
  if (log(stats::runif(1, 0, 1)) < alpha1) {
    k <- k.star  # Accept k.star
    f.store <- f.k.star
  }
  #####
  # End: Step 1
  #####
  
  #####
  # Step 2: Metropolis-within-Gibbs step for V (EXPENSIVE)
  #####
  for (l in 1:L) {
    
    V.star <- V.old <- V
    if (l > 1) {
      for (il in 1:(l - 1)) {
        V.star[il] <- V.old[il] <- V[il]
      }
    }
    
    # Uniform proposal (V[,i] - eps, V[,i] + eps) on (0,1)
    V.star[l] <- stats::runif(1, V.star[l] - eps[l], V.star[l] + eps[l])
    V.star[l][V.star[l] > 1] <- V.star[l] - 1  # Puts in [0, 1]
    V.star[l][V.star[l] < 0] <- V.star[l] + 1  # Puts in [0, 1]
    
    # log posterior for proposal
    f.V.star <- lpost(omega, FZ, V.star, W, k, tau,
                      M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                      pdgrm, db.list)   

    # log posterior of previous iteration
    f.V <- f.store
    
    # Accept/reject
    alpha2 <- min(0, f.V.star - f.V)  # log acceptance ratio
    if (log(stats::runif(1, 0, 1)) < alpha2) {
      V[l] <- V.star[l]  # Accept V.star
      f.store <- f.V.star
    }
    
  }  
  #####
  # End: Step 2
  #####
  
  #####
  # Step 3: Metropolis-within-Gibbs step for W (EXPENSIVE)
  #####
  for (l in 1:(L + 1)) {
    
    W.star <- W.old <- W
    if (l > 1) {
      for (il in 1:(l - 1)) {
        W.star[il] <- W.old[il] <- W[il]
      }
    }
    
    # Uniform proposal from (W[,i] - eps, W[,i] + eps) on (0,1)
    W.star[l] <- stats::runif(1, W.star[l] - eps[l], W.star[l] + eps[l])
    W.star[l][W.star[l] > 1] <- W.star[l] - 1  # Puts in [0, 1]
    W.star[l][W.star[l] < 0] <- W.star[l] + 1  # Puts in [0, 1]
    
    # log posterior for proposal
    f.W.star <- lpost(omega, FZ, V, W.star, k, tau,
                      M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                      pdgrm, db.list)   
    
    # log posterior for previous iteration
    f.W <- f.store
    
    # Accept/reject
    alpha3 <- min(0, f.W.star - f.W)  # log acceptance ratio
    if(log(stats::runif(1, 0, 1)) < alpha3) {
      W[l] <- W.star[l]  # Accept W.star
      f.store <- f.W.star
    }
    
  }  
  #####
  # End: Step 3
  #####
  
  #####
  # Step 4: Directly sample tau from conjugate Inverse-Gamma density
  #####
  q.psd <- qpsd(omega, V, W, k, db.list)$psd
  q <- unrollPsd(q.psd, n)
  
  # Note: (n - 1) and (n - 2) here.  Remove the first and last terms for even and first for odd
  if (n %% 2) {  # Odd length series
    bFreq <- 1
    tau <- 1 / stats::rgamma(1, tau.alpha + (n - 1) / 2, 
                             tau.beta + sum(pdgrm[-bFreq] / q[-bFreq]) / (2 * pi) / 2)
  }
  else {  # Even length series
    bFreq <- c(1, n)
    tau <- 1 / stats::rgamma(1, tau.alpha + (n - 2) / 2, 
                             tau.beta + sum(pdgrm[-bFreq] / q[-bFreq]) / (2 * pi) / 2)
  }
  #####
  # End: Step 4
  #####
  
  # Output list of parameters sampled
  return(list(k = k,
              V = V,
              W = W,
              tau = tau,
              q = q))
  
}
