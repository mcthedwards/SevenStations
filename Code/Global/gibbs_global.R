# Sample the global slope (and intercept) parameters
gibbs_global = function(beta, beta.v0, mu.0, sigma2.0) {
  
  s = length(beta)  # Should be 7 for 7 Stations 
                    # Works otherwise for 1 Station
  
  sigma2.n = 1 / (s / beta.v0 + 1 / sigma2.0)
  mu.n = (sum(beta) / beta.v0 + mu.0 / sigma2.0) * sigma2.n
  beta.G = rnorm(1, mu.n, sqrt(sigma2.n))  # Sample global slope (and intercept)
  
  return(beta.G)
  
}