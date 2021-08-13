# Bayesian estimation of the long-term linear trend of New Zealand annual average temperatures

An accurate understanding of the long-term evolution of the temperature is key to understanding the impact of global warming. NIWA provides time series of annual average temperatures at various different sites in New Zealand.
However, these temperature time series lack homogeneity due to changes in instrumentation and re-siting of recording stations that has necessitated adjustments in the past. The goal of this project is to use a hierarchical Bayesian model to estimate the slope. It will aim for a robust analysis by using a nonparametric approach to model the time series errors. The first phase of the project will be concerned with data wrangling, accessing the data from the NIWA website and bringing it into a suitable format for subsequent analysis using R packages for nonparametric time series errors. Exploratory analysis will be part of this phase. In a second phase, the model for time series errors will need to be combined with a hierarchical linear model for the slope. Sampling from the posterior distribution will either be performed using JAGS and/or Metropolis-Hastings routines written in R.

Requirements: A good knowledge of and interest in Bayesian inference, MCMC techniques, and time series as well as good programming skills and knowledge of R and JAGS are essential.

## Data
Can be downloaded directly from CliFlo (https://cliflo.niwa.co.nz/).  You can also use the clifro R package.  

## Model
We can think of the model as Data = Signal + Noise.  We will use a blocked Gibbs sampler to sample the signal parameters given the noise parameters, and the noise parameters given the signal parameters.   

## Noise model
Whittle likelihood with Bernstein-Dirichlet prior (Choudhuri et al. 2004).  Could also consider the corrected likelihood with Bernstein-Dirichlet prior (Kirch et al. 2019), the B-spline prior (Edwards et al. 2019), or P-spline prior (Maturana-Russel and Meyer 2021).

## Signal model
Linear regression with time as the explanatory variable, hierarchical slope for annual temperature at the seven stations, different y-intercepts to account for (known) level-shifts.

## Progress:
Implemented a 2-station (Auckland and Wellington) version of the model, with same slope, and different intercepts for each (known) weather station change-point/level shift.

## Code:
- We require the following R libraries: Rcpp, compiler, MASS.
- gibbs_functions.R: All the functions used in the background, e.g., FFT, log likelihood, log posterior, stick-breaking, plot methods, etc.  
- gibbs_noise.R: The blocked noise sampler.  Outputs an iteration of noise parameters.
- gibbs_sampler_noise_only.R: The MCMC sampler if we are only interested in estimating the PSD of a noise time series.  Function is gibbs_psd().
- gibbs_signal.R: The blocked signal sampler.  It outputs an iteration of regression coefficients.  Note: May want to implement Cholesky version of this.
- gibbs_sampler_with_signal.R: The blocked Gibbs MCMC sampler that iterates between noise parameters and signal parameters.  Function is gibbs_temperature().  This is the code we need to change to implement a hierarchical version of the model.
- gibbs_example.R: Two examples. First one is a noise-only run on simulated AR data. Second is a blocked sampler where we have Auckland annual average temperature data from 1910 onwards.  We want to estimate a linear slope.  There are four levels due to the station changing over time.  This is represented by using four different intercepts in the design matrix.  The noise parameters are also estimated and we can construct the PSD if we want, although we are primarily interested in the signal.

## To do:
- Tidy up old code - DONE
- Scrape most recent data from NIWA - DONE
- Implement hierarchical version for linear slope coefficient
- Impute missing data. Originally used an ad hoc fill-in for the one missing Wellington data point. Implemented in beyondWhittle R package for the noise component. 

