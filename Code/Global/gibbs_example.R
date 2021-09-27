source("gibbs_functions.R")
source("gibbs_noise.R")
source("gibbs_signal.R")
source("gibbs_global.R")
source("gibbs_sampler_with_signal.R")

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
                         20000, 10000, 1)

# Plot signal estimate
plot(mcmc, legend.loc = NA)

plot.ts(mcmc$beta)
plot.ts(mcmc$beta.G)
plot.ts(mcmc$alpha)
