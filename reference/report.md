# Get reported quantities from and `RTMB` object and return a `LaMaModel`

Having fitted a latent Markov model using automatic differentiation via
`RTMB`, this function calls `RTMB`'s `obj$report()` and does some
additional processing. This then yields estimated parameters on their
natural scale, allows for convenient calculation of `AIC` and `BIC`,
state-decoding, and residual calculation.

## Usage

``` r
report(obj)
```

## Arguments

- obj:

  Object returned by `RTMB::MakeADFun()`

## Value

A model object of class "`LaMaModel`" containing a list with the
reported quantities from the `RTMB` object, along with log-likelihood,
number of parameters, and number of observations.

## Examples

``` r
data <- trex[1:200,]

# initial parameters and observations
par <- list(
  log_mu = log(c(0.3, 1)),      # initial means for step length (log-transformed)
  log_sigma = log(c(0.2, 0.7)), # initial sds for step length (log-transformed)
  eta = rep(-2, 2)              # initial t.p.m. parameters (on logit scale)
)    
dat <- list(
  step = data$step,             # hourly step lengths
  nStates = 2                   # number of hidden states
)

# likelihood function
nll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta)
  delta <- stationary(Gamma)
  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  allprobs <- matrix(1, length(step), nStates)
  ind <- which(!is.na(step))
  for(j in 1:nStates) {
    allprobs[ind,j] <- dgamma2(step[ind], mu[j], sigma[j]) 
  } 
  -forward(delta, Gamma, allprobs)
}

# automatic differentiation and optimisation
obj <- MakeADFun(nll, par, silent = TRUE)
#> Performance tip: Consider running `TapeConfig(matmul = 'plain')` before `MakeADFun()` to speed up the forward algorithm.
opt <- nlminb(obj$par, obj$fn, obj$gr)

### reporting ###
mod <- report(obj)

# estimated quantities on natural scale
mod$mu
#> [1] 0.3136312 2.5678297
mod$sigma
#> [1] 0.2033994 1.5138886
mod$Gamma
#>           S1        S2
#> S1 0.8085298 0.1914702
#> S2 0.2255601 0.7744399

# information criteria
AIC(mod)
#> [1] 432.5328
BIC(mod)
#> [1] 452.3227

# state decoding
states <- viterbi(mod = mod)   # global decoding
probs <- stateprobs(mod = mod) # local decoding

# residual calculation
pres <- pseudo_res(data$step, # observation sequence
                   "gamma2",  # distribution family
                   list(mean = mod$mu, sd = mod$sigma), # parameters for that family
                   mod = mod) # model object
```
