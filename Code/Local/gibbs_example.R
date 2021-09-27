source("gibbs_functions.R")
source("gibbs_noise.R")
source("gibbs_sampler_noise_only.R")
source("gibbs_signal.R")
source("gibbs_sampler_with_signal.R")

#####
# 1. Noise-only MCMC
#####

# Simulate data
data = arima.sim(n = 512, model = list(ar = c(0.9, -0.9)))

# Run MCMC
mcmc = gibbs_psd(data, 10000, 5000, 1)

# Plot spectrum estimate
plot(mcmc)


#####
# 2. Signal plus noise MCMC
#####

# Read in Auckland temperature data
data = read.csv("annual_temp_2016.csv")
data = data[data$location=="Auckland", ]

times = data$year
temp = data$temp

# Change points
cp = c()
cp[1] = 1
for (i in 2:nrow(data)) {
  cp[i] = (data$site_number[i] != data$site_number[i - 1])
}
cp = which(cp == 1)  # Auckland change points

# Run MCMC
mcmc = gibbs_temperature(temp, times, cp,
                         10000, 5000, 1)

# Plot signal estimate
plot(mcmc, legend.loc = NA)


