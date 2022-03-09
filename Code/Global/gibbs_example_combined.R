source("gibbs_functions.R")
source("gibbs_noise.R")
source("gibbs_signal.R")
source("gibbs_global.R")
source("gibbs_combined.R")


#old data
#data = read.csv("annual_temp_2016.csv")

data = read.csv("updated_data_2021.csv")

sapply(data,class)
data[is.na(data$temp),]

#replace NAs with previous non-null value
data$temp[is.na(data$temp)] = data$temp[which(is.na(data$temp))-1]
data[is.na(data$temp),]


#three NA's in Lincoln 1916-1918 - replace last two with 12. Aso Lincoln in 1992 - replace with 12
data$temp[is.na(data$temp)] = 12

data[is.na(data$temp),]

#try without Lincoln
six_stations_data = data[data$location != "Lincoln",]

# Run MCMC
mcmc = gibbs_temperature_overall(six_stations_data, 20000, 10000, nloc = 6)

# Plot signal estimate
# plot(mcmc, legend.loc = NA)

plot.ts(mcmc$beta.G)
plot.ts(mcmc$betas[[1]])

mean(mcmc$beta.G)
mean(mcmc$beta.G)*100
sapply(1:6, function (x) mean(mcmc$betas[[x]]))*100

sapply(1:6, function (x) quantile(mcmc$betas[[x]], probs = c(0.05, 0.5, 0.95)))*100
quantile(mcmc$beta.G, probs = c(0.05, 0.5, 0.95)) * 100

#global slope samples still distributed around zero. Need to fix. 

