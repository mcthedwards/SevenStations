# Bayesian estimation of the long-term linear trend of New Zealand annual average temperatures

An accurate understanding of the long-term evolution of the temperature is key to understanding the impact of global warming. NIWA provides time series of annual average temperatures at various different sites in New Zealand.
However, these temperature time series lack homogeneity due to changes in instrumentation and re-siting of recording stations that has necessitated adjustments in the past. The goal of this project is to use a hierarchical Bayesian model to estimate the slope. It will aim for a robust analysis by using a nonparametric approach to model the time series errors. The first phase of the project will be concerned with data wrangling, accessing the data from the NIWA website and bringing it into a suitable format for subsequent analysis using R packages for nonparametric time series errors. Exploratory analysis will be part of this phase. In a second phase, the model for time series errors will need to be combined with a hierarchical linear model for the slope. Sampling from the posterior distribution will either be performed using JAGS and/or Metropolis-Hastings routines written in R.

Requirements: A good knowledge of and interest in Bayesian inference, MCMC techniques, and time series as well as good programming skills and knowledge of R and JAGS are essential.

## Noise model
Corrected likelihood with Bernstein-Dirichlet prior (Kirch et al. 2019).

## Signal model
Linear regression with time as the explanatory variable, hierarchical slope for annual temperature at the seven stations, different y-intercepts to account for (known) level-shifts.

## Progress:
Implemented a 2-station (Auckland and Wellington) version of the model, with same slope, and different intercepts for each (known) weather station change-point/level shift.

## To do:
- Find most recent version of the code
- Scrape most recent data from NIWA
- Implement hierarchical version for linear slope coefficient
- Explore P-spline or B-spline noise model instead? Perhaps simplify to original Bernstein-Dirichlet prior with Whittle likelihood.
- Figure out how to impute missing data. Originally used an ad hoc fill-in for the one missing Wellington data point.

