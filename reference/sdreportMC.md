# Monte Carlo version of `sdreport`

After optimisation of an AD model, `sdreportMC` can be used to calculate
samples of confidence intervals of all model parameters and
[`REPORT()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)ed
quantities including nonlinear functions of random effects and
parameters.

## Usage

``` r
sdreportMC(
  obj,
  what,
  nSamples = 1000,
  Hessian = NULL,
  CI = FALSE,
  probs = c(0.025, 0.975)
)
```

## Arguments

- obj:

  object returned by `MakeADFun()` after optimisation or model of class
  `qremlModel` as returned by
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md).

- what:

  vector of strings with names of parameters and/or
  [`REPORT()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)ed
  quantities to be reported

- nSamples:

  number of samples to draw from the multivariate normal distribution of
  the MLE

- Hessian:

  optional Hessian matrix. If not provided, it will be computed from the
  object

- CI:

  logical. If `TRUE`, only confidence intervals instead of samples will
  be returned

- probs:

  vector of probabilities for the confidence intervals (ignored if no
  CIs are computed)

## Value

named list corresponding to the elements of `what`. Each element has the
structure of the corresponding quantity with an additional dimension
added for the samples. For example, if a quantity is a vector, the list
contains a matrix. If a quantity is a matrix, the list contains an
array. If quantity is an array, the list contains an array with one
extra dimension.

## Details

This function simply samples from the approximate multivariate normal
distribution of the maximum likelihood estimate (MLE) of the parameters
\$\$ \hat{\theta} \sim N(\theta, H^{-1}), \$\$ where \\H\\ is the
Hessian matrix of the negative log-likelihood function at the MLE. It
then returns either the sampled parameters or
[`REPORT()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html)ed
transformations of them. If `CI = TRUE`, it does not return the samples
but directly computes confidence intervals.

If you are interested in several quantities, calling `sdreportMC` once
with a vector `what` will generally be faster than calling it several
times with single elements of `what`.

## Examples

``` r
# fitting an HMM to the trex data and running sdreportMC
## negative log-likelihood function
nll = function(par) {
  getAll(par, dat) # makes everything contained available without $
  Gamma = tpm(eta) # computes transition probability matrix from unconstrained eta
  delta = stationary(Gamma) # computes stationary distribution
  # exponentiating because all parameters strictly positive
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  # reporting statements for sdreportMC
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle)) # only for non-NA obs.
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])*dvm(angle[ind],0,kappa[j])
  }
  -forward(delta, Gamma, allprobs) # simple forward algorithm
}

## initial parameter list
par = list(
 logmu = log(c(0.3, 1)),       # initial means for step length (log-transformed)
  logsigma = log(c(0.2, 0.7)), # initial sds for step length (log-transformed)
  logkappa = log(c(0.2, 0.7)), # initial concentration for turning angle (log-transformed)
  eta = rep(-2, 2)             # initial t.p.m. parameters (on logit scale)
)   
## data and hyperparameters
dat = list(
  step = trex$step[1:500],   # hourly step lengths
  angle = trex$angle[1:500], # hourly turning angles
  N = 2
)

## creating AD function
obj = MakeADFun(nll, par, silent = TRUE) # creating the objective function

## optimising
opt = nlminb(obj$par, obj$fn, obj$gr) # optimization

## running sdreportMC
# `mu` has report statement, `delta` is automatically reported by `forward()`
sdrMC = sdreportMC(obj, 
                   what = c("mu", "delta"), 
                   nSamples = 50)
#> Sampling reported quantities...
dim(sdrMC$delta)
#> [1] 50  2
# now a matrix with 50 samples (rows)
```
