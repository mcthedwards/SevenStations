Noise model:  Corrected likelihood with Bernstein-Dirichlet prior (Kirch et al. 2019).

Signal model:  Linear regression with time as the explanatory variable, hierarchical slope for annual temperature at the seven stations, different y-intercepts to account for (known) level-shifts.

Done:
Implemented a 2-station (Auckland and Wellington) version of the model, with same slope, and different intercepts for each (known) weather station change-point/level shift.

To do:
- Find most recent version of the code
- Scrape most recent data from NIWA
- Implement hierarchical version for linear slope coefficient
- Explore P-spline or B-spline noise model instead? Perhaps simplify to original Bernstein-Dirichlet prior with Whittle likelihood.
- Figure out how to impute missing data. Originally used an ad hoc fill-in for the one missing Wellington data point.

