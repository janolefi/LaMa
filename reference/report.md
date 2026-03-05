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
data <- trex[1:500,]

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
  REPORT(step)
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
obj <- MakeADFun(nll, par)
#> Performance tip: Consider running `TapeConfig(matmul = 'plain')` before `MakeADFun()` to speed up the forward algorithm.
opt <- nlminb(obj$par, obj$fn, obj$gr)
#> outer mgc:  557.295 
#> outer mgc:  351.0983 
#> outer mgc:  287.616 
#> outer mgc:  47.5427 
#> outer mgc:  14.6837 
#> outer mgc:  27.66109 
#> outer mgc:  9.562208 
#> outer mgc:  9.821129 
#> outer mgc:  13.44982 
#> outer mgc:  19.303 
#> outer mgc:  35.85799 
#> outer mgc:  7.582501 
#> outer mgc:  4.230678 
#> outer mgc:  4.698342 
#> outer mgc:  5.808632 
#> outer mgc:  3.698382 
#> outer mgc:  4.58231 
#> outer mgc:  7.42501 
#> outer mgc:  10.8721 
#> outer mgc:  1.829161 
#> outer mgc:  4.166831 
#> outer mgc:  1.591044 
#> outer mgc:  2.216416 
#> outer mgc:  1.953366 
#> outer mgc:  0.4959815 
#> outer mgc:  1.196605 
#> outer mgc:  0.3334328 
#> outer mgc:  0.3744078 
#> outer mgc:  0.3096951 
#> outer mgc:  0.4332586 
#> outer mgc:  0.08363097 
#> outer mgc:  0.007251231 
#> outer mgc:  0.0002769499 
#> outer mgc:  1.449259e-05 

### reporting ###
mod <- report(obj)

# estimated quantities on natural scale
mod$mu
#> [1] 0.3288388 2.4972538
mod$sigma
#> [1] 0.2275702 1.3900984
mod$Gamma
#>           S1        S2
#> S1 0.8284923 0.1715077
#> S2 0.1789212 0.8210788

# information criteria
AIC(mod)
#> [1] 1084.425
BIC(mod)
#> [1] 1109.713

# state decoding
states <- viterbi(mod = mod)   # global decoding
probs <- stateprobs(mod = mod) # local decoding

# residual calculation
pres <- pseudo_res(data$step, # observation sequence
                   "gamma2",  # distribution family
                   list(mean = mu, sd = sigma), # parameters for that family
                   mod = mod) # model object
#> Error: object 'mu' not found
```
