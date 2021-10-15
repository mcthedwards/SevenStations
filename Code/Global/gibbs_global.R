# Sample the global slope parameters
gibbs_global_slope = function(beta, beta.v0, mu.0, sigma2.0) {
  
  # Normal conjugate (hyper)prior
  
  s = length(beta)  # Should be 7 for 7 Stations 
                    # Works otherwise for 1 Station
  
  sigma2.n = 1 / (s / beta.v0 + 1 / sigma2.0)
  mu.n = (sum(beta) / beta.v0 + mu.0 / sigma2.0) * sigma2.n
  beta.G = rnorm(1, mu.n, sqrt(sigma2.n))  # Sample global slope 
  
  return(beta.G)
  
}

# Sample global variance
gibbs_global_variance = function(beta, beta.G, a0, b0) {
  
  # Inverse-Gamma conjugate (hyper)prior
  
  s = length(beta)  # Should be 7 for 7 Stations 
  
  a.n = a0 + s / 2
  b.n = b0 + 0.5 * sum((beta - beta.G) ^ 2)
  beta.v0 = 1 / rgamma(1, a.n, b.n)  # Sample global variance
  
  return(beta.v0)
  
}