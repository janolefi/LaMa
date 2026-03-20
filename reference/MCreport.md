# Sample parameters from approximate Gaussian posterior distribution

Efficient Monte Carlo sampling of parameters (and `REPORT`ed quantities)
from the approximate posterior of an `RTMB` model. See
[sdreport](https://rdrr.io/pkg/TMB/man/sdreport.html) for details on
posterior variance-covariance in random effects models.

## Usage

``` r
MCreport(obj, nSamples = 1000, include_random_pars = TRUE, report = FALSE, ...)
```

## Arguments

- obj:

  Optimised `RTMB` object

- nSamples:

  Number of samples to draw

- include_random_pars:

  Logical; Should random parameters be included in the output?

- report:

  Logical; Should `REPORT`ed quantities be sampled as well? Defaults to
  `FALSE` because this may be slow depending on your model.

- ...:

  For internal use only

## Value

A list structured like the original parameter list used in the MakeADFun
call (potentially including additional `REPORT`ed quantities). Each
entry is a list with `nSamples` entries.

## Examples

``` r
step <- trex$step[1:1000] # subsetting trex data
N <- 2                    # 2 states

# custom likelihood
nll <- function(par) {
  getAll(par)
  Gamma <- tpm(eta)
  delta <- stationary(Gamma)
  mu <- exp(log_mu); REPORT(mu)
  sigma <- exp(log_sigma); REPORT(sigma)
  allprobs <- matrix(1, length(step), N)
  for(j in 1:N) allprobs[,j] <- dgamma2(step, mu[j], sigma[j])
  -forward(delta, Gamma, allprobs)
}

# initial parameters in named list
par0 <- list(eta = rep(-2,2), 
             log_mu = log(c(0.3, 1)), 
             log_sigma = log(c(0.2, 0.7)))
    
# constructing AD object        
obj <- MakeADFun(nll, par0, silent = TRUE)
#> Performance tip: Consider running TapeConfig(matmul = 'plain') before MakeADFun() to speed up the forward algorithm.

# optimising
opt <- nlminb(obj$par, obj$fn, obj$gr)

# sampling from distribution of the MLE
samples <- MCreport(obj, nSamples = 10, report = TRUE)
#> Evaluating Hessian...
#> Computing reported quantities...
```
