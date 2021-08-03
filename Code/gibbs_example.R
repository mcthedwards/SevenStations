source("gibbs_functions.R")
source("gibbs_noise.R")
source("gibbs_run.R")

# Simulate data
data = arima.sim(n = 256, model = list(ar = 0.9))

# Run MCMC
mcmc = gibbs_run(data, 10000, 5000, 1)

# Plot spectrum estimate
plot(mcmc)


