library(Rcpp)
library(compiler)

# C++ function for generating  p from v in Stick Breaking DP representation
cppFunction("NumericVector pFromV(NumericVector v) {
  unsigned L = v.size();
  NumericVector p(L + 1);
  double currentProduct = 1.0;
  double pSum = 0.0;
  for (unsigned l = 0; l < L; ++l) {
    p[l + 1] = currentProduct * v[l];
    currentProduct *= (1.0 - v[l]);
    pSum += p[l + 1];
  }
  p[0] = std::max(1.0 - pSum, 0.0); // account for numerical instabilities
  return p;
}")

# C++ function for generating  v from p (inverse stick breaking)
# NOTE: p is assumed to have length L, i.e. it does NOT contain p_0 !!
cppFunction("NumericVector vFromP(NumericVector p, const double eps=1e-8) {
  unsigned L = p.size();
  NumericVector v(L);
  double currentProduct = 1.0;
  for (unsigned l = 0; l < L; ++l) {
    v[l] = std::min(std::max(p[l] / currentProduct, eps),1.0-eps); // numerical stability
    //v[l] = p[l] / currentProduct;
    currentProduct *= (1.0 - v[l]);
  }
  return v;
}")

# C++ function for computing mixture weights of Bernstein-Mixtures given the probabilities p, values w, and degree k
cppFunction("NumericVector mixtureWeight(NumericVector p, NumericVector w, unsigned k) {
  typedef std::pair<double, double> wpType;
  std::vector<wpType> wp;
  for (unsigned l = 0; l < p.size(); ++l) {
    wp.push_back(wpType(w[l], p[l]));
  }
  std::sort(wp.begin(), wp.end());
  NumericVector weight(k);
  unsigned l = 0;
  for (unsigned j = 1; j <= k; ++j) {
    weight[j-1] = 0;
    double wMax = j / (double)k;
    while (l < wp.size() && wp[l].first <= wMax) {
      weight[j-1] += wp[l].second;
      l += 1;
    }
  }
  return weight;
}")

# C++ function for building a density mixture, given mixture weights and functions
cppFunction("NumericVector densityMixture(NumericVector weights, NumericMatrix densities) {
  if (weights.size() != densities.nrow()) {
    return(NumericVector());
  }
  const unsigned n = densities.ncol();
  NumericVector res(n);
  for (unsigned omega = 0; omega < n; ++omega) {
    res[omega] = 0.0;
  }
  for (unsigned j = 0; j < weights.size(); ++j) {
    for (unsigned omega = 0; omega < n; ++omega) {
      res[omega] += weights[j] * densities(j, omega);
    }
  }
  return(res);
}")

# C++ help function to redundantly roll out a PSD to length n
cppFunction("NumericVector unrollPsd(NumericVector qPsd, unsigned n) {
  NumericVector q(n);
  q[0] = qPsd[0];
  const unsigned N = (n-1)/2;
  for (unsigned i = 1; i <= N; ++i) {
    const unsigned j = 2 * i - 1;
    q[j] = qPsd[i];
    q[j+1] = qPsd[i];
  }
  if (!(n % 2)) {
    q[n-1] = qPsd[qPsd.size() - 1];
  }
  return(q);
}")

# Compute a PSD in the Bernstein-Dirichlet parametrization
qpsd <- function(omega, v, w, k, db.list, epsilon=1e-20) {
  p <- pFromV(v)
  weight <- mixtureWeight(p, w, k)
  psd <- densityMixture(weight, db.list[[k]])
  psd <- pmax(psd, epsilon)
  return(list(psd = psd,
              weight = weight,
              p=p)) 
}

# log Whittle likelihood
llike <- function(omega, FZ, v, w, k, tau, pdgrm, db.list) {
  
  # Calculates Whittle log-likelihood for Gaussian errors
  n <- length(FZ)
  
  # Un-normalised PSD (defined on [0, 1])
  q.psd <- qpsd(omega, v, w, k, db.list)$psd
  q <- unrollPsd(q.psd, n)
  
  # Normalised PSD (defined on [0, pi])
  f <- tau * q
  
  # Do not consider the following frequencies in the likelihood
  # They correspond to the mean (or alternating mean)
  if (n %% 2) {  # Odd
    bFreq <- 1
  } else {  # Even
    bFreq <- c(1, n)
  }
  
  # Whittle log-likelihood
  llike <- -0.5 * sum(log(f[-bFreq] * 2 * pi) + pdgrm[-bFreq] / (f[-bFreq] * 2 * pi)) 
  
  return(llike)
  
}

# log prior of Bernstein-Dirichlet mixture - all unnormalized
lprior <- function(v, w, k, tau, M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta) {
  
  # log joint prior - all unnormalised
  # Hyperparameters are M, g0.a, g0.b, k.theta, tau.alpha, tau.beta

  lp <- (M - 1) * sum(log(1 - v)) +  # log prior for V's - beta(1, M)
    sum((g0.alpha - 1) * log(w) + (g0.beta - 1) * log(1 - w)) -  # log prior for Z's - beta(a, b)
    k.theta * k * log(k) -   # log prior for k
    (tau.alpha + 1) * log(tau) - tau.beta / tau # log prior for tau (Inverse Gamma)
  
  return(lp)
  
}

# log posterior
lpost <- function(omega, FZ, v, w, k, tau,
                  M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta,
                  pdgrm, db.list) {
  
  # Unnormalised log posterior
  ll <- llike(omega, FZ, v, w, k, tau, pdgrm, db.list)
  lp <- lprior(v, w, k, tau, M, g0.alpha, g0.beta, k.theta, tau.alpha, tau.beta)
  
  return(ll+lp)
  
}

# Compute F_n X_n with the real-valued Fourier matrix F_n
# Function computes FZ (i.e. fast Fourier transformed data)
# Outputs coefficients in correct order and rescaled
fast_ft <- compiler::cmpfun(function(x) {
  
  n <- length(x)
  sqrt2 <- sqrt(2)
  sqrtn <- sqrt(n)
  
  # Cyclically shift so last observation becomes first
  x <- c(x[n], x[-n])  # Important since fft() uses 0:(n-1) but we use 1:n
  
  # FFT
  fourier <- fft(x)
  
  # Extract non-redundant real and imaginary coefficients in correct order and rescale
  FZ <- rep(NA, n)
  FZ[1] <- Re(fourier[1]) # First coefficient
  if (n %% 2) {  # Even
    N <- (n - 1) / 2
    FZ[2 * (1:N)] <- sqrt2 * Re(fourier[2:(N + 1)]) # Real coefficients
    FZ[2 * (1:N) + 1] <- sqrt2 * Im(fourier[2:(N + 1)]) # Imaginary coefficients
  } else {  # Odd
    FZ[n] <- Re(fourier[n / 2 + 1]) # Last coefficient
    FZ[2 * 1:(n / 2 - 1)] <- sqrt2 * Re(fourier[2:(n / 2)]) # Real coefficients
    FZ[2 * 1:(n / 2 - 1) + 1] <- sqrt2 * Im(fourier[2:(n / 2)]) # Imaginary coefficients
  }
  
  return(FZ / sqrtn)
  
})

plot.psd = function(x, legend.loc = "topright", ylog = TRUE, ...) {  # Plot method for "psd" class
  N = length(x$pdgrm)
  freq = seq(0, pi, length = N)
  # Frequencies to remove from estimate
  if (x$n %% 2) {  # Odd length time series
    bFreq = 1
  }
  else {  # Even length time series
    bFreq = c(1, N)
  }
  
  if (!is.logical(ylog)) stop("ylog must be TRUE or FALSE")
  
  if (ylog == TRUE) {
    graphics::plot.default(freq[-bFreq], log(x$pdgrm[-bFreq]), type = "l", col = "grey",
                           xlab = "Frequency", ylab = "log PSD", ...)
    graphics::lines(freq[-bFreq], log(x$psd.median[-bFreq]), lwd = 2)
    graphics::lines(freq[-bFreq], log(x$psd.p05[-bFreq]), lwd = 2, lty = 2, col = 4)
    graphics::lines(freq[-bFreq], log(x$psd.p95[-bFreq]), lwd = 2, lty = 2, col = 4)
    if (!is.na(legend.loc)) {
      graphics::legend(legend.loc, legend = c("periodogram", "posterior median", "90% credible region"), 
                       col = c("grey", "black", "blue"), lwd = c(1, 2, 2), lty = c(1, 1, 2))  
    }
  }
  if (ylog == FALSE) {
    graphics::plot.default(freq[-bFreq], x$pdgrm[-bFreq], type = "l", col = "grey",
                           xlab = "Frequency", ylab = "PSD", ...)
    graphics::lines(freq[-bFreq], x$psd.median[-bFreq], lwd = 2)
    graphics::lines(freq[-bFreq], x$psd.p05[-bFreq], lwd = 2, lty = 2, col = 4)
    graphics::lines(freq[-bFreq], x$psd.p95[-bFreq], lwd = 2, lty = 2, col = 4)
    if (!is.na(legend.loc)) {
      graphics::legend(legend.loc, legend = c("periodogram", "posterior median", "90% credible region"), 
                       col = c("grey", "black", "blue"), lwd = c(1, 2, 2), lty = c(1, 1, 2))  
    }
  }
}