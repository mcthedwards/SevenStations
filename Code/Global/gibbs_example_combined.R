source("gibbs_functions.R")
source("gibbs_noise.R")
source("gibbs_signal_global.R")
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


# Run MCMC
mcmc = gibbs_temperature_overall(data, 20000, 10000)

# Plot signal estimate
# plot(mcmc, legend.loc = NA)

plot.ts(mcmc$beta.G)
plot.ts(mcmc$betas[[1]])
plot.ts(mcmc$alphas[[1]])

mean(mcmc$beta.G)
mean(mcmc$beta.G)*100
sapply(1:7, function (x) mean(mcmc$betas[[x]]))*100
hist(mcmc$beta.G)


#global slope samples still distributed around zero. Need to fix. 

