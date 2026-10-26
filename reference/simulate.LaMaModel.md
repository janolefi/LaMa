# Simulate from a fitted model

Simulates states and observations from a model fitted with
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md) and
[`statedens`](https://janolefi.github.io/LaMa/reference/statedens.md).
By default, this is done on the observed design: the tracks and their
lengths, the initial distribution `delta`, the transition probabilities
`Gamma` and the parameters of the state-dependent densities are those
reported by the model, i.e. conditional on the observed covariates.
Simulating for new covariates is not possible; instead, rebuild the
likelihood function for the new data and call
[`lama`](https://janolefi.github.io/LaMa/reference/lama.md) on the
resulting object.

## Usage

``` r
# S3 method for class 'LaMaModel'
simulate(
  object,
  nsim = 1,
  seed = NULL,
  rdens = NULL,
  random = FALSE,
  fill_na = FALSE,
  ...
)
```

## Arguments

- object:

  model object from
  [`lama`](https://janolefi.github.io/LaMa/reference/lama.md) or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) that was
  built with
  [`statedens`](https://janolefi.github.io/LaMa/reference/statedens.md)

- nsim:

  number of simulations

- seed:

  optional seed for the random number generator

- rdens:

  optional named list of random number generators of the data streams,
  for densities without one named by replacing the leading `d` with `r`
  (e.g. from
  [`zero_inflate`](https://janolefi.github.io/LaMa/reference/zero_inflate.md)).
  The first argument is the number of draws, the others are the
  parameters as in the density.

- random:

  logical; if `TRUE`, `obj$simulate()` of `RTMB` is used, so that random
  effects are drawn from their distribution first, see ‘Details’. The
  data streams are those passed to
  [`statedens`](https://janolefi.github.io/LaMa/reference/statedens.md).
  Default is `FALSE`, which keeps the reported (estimated) random
  effects.

- fill_na:

  logical; if `TRUE`, missing observations of the data are filled in by
  simulated values (their parameters need to be defined, e.g. covariates
  must not be missing). Default is `FALSE`, which keeps them missing.
  Streams that are missing completely are always simulated.

- ...:

  currently not used

## Value

data frame with the columns `trackID` (index of the track, if used),
`state` and one column per data stream, for `nsim = 1`, else a list of
`nsim` such data frames.

## Details

The states are drawn as in the decoding functions, so that all types of
`Gamma` that
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md)
accepts (homogeneous, track-specific, time-varying, periodic) work, as
does a continuous-time model with one transition matrix per time
interval (the simulation is then conditional on the observed time
points). Each data stream is drawn from the random number generator
belonging to its density, e.g. `rgamma2` for `dgamma2`.

Missing observations of the data stay missing (unless `fill_na = TRUE`),
except for a stream that is missing completely: such a placeholder (e.g.
`step = rep(NA, n)`) lets you simulate from a likelihood function
without data, see the examples.

Models with random effects are by default simulated conditional on the
estimated random effects. With `random = TRUE`, `RTMB` draws the random
effects first, which requires that the density of each random effect is
evaluated in the likelihood function *before* the random effect is used,
either by `u %~% dnorm(...)` or by
`nll <- nll - sum(dnorm(u, ..., log = TRUE))` (or another density). The
states and data are then simulated given these draws, so that every
simulation has new random effects. This also requires that `allprobs` is
the product of the data streams built by
[`statedens`](https://janolefi.github.io/LaMa/reference/statedens.md),
with each data stream passed once.

## Examples

``` r
nll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta)
  allprobs <- statedens(step, dgamma2, exp(log_mu), exp(log_sigma)) *
    statedens(angle, dvm, 0, exp(log_kappa))
  -forward(stationary(Gamma), Gamma, allprobs)
}
dat <- list(step = trex$step[1:500], angle = trex$angle[1:500])
par <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)),
            log_kappa = log(c(0.2, 1.5)))
obj <- MakeADFun(nll, par, silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
mod <- lama(obj)

sim <- simulate(mod, seed = 1)
head(sim)
#>   state      step      angle
#> 1     1 0.2800839         NA
#> 2     1 0.3407063 -1.0075237
#> 3     1 0.2379878  0.6967544
#> 4     2 1.3826001 -0.4879542
#> 5     2 1.2731907 -3.0053890
#> 6     1 0.5038482 -1.4462509

# simulate without data: all-NA streams as placeholders, parameters of your choice, then refit
dat$step <- dat$angle <- rep(NA, 500)
obj0 <- MakeADFun(nll, par, silent = TRUE)
sim <- simulate(lama(obj0), seed = 1)
dat$step <- sim$step; dat$angle <- sim$angle
obj1 <- MakeADFun(nll, par, silent = TRUE)
opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)

# random effects: RTMB draws them first, then the states and data are simulated given them
nll_re <- function(par) {
  getAll(par, dat)
  u %~% dnorm(0, exp(log_su)) # track-level random effect on the mean step length in state 1
  mu <- exp(cbind(log_mu[1] + u[id], log_mu[2]))
  Gamma <- tpm(eta)
  .nll <- .nll - forward(stationary(Gamma), Gamma, statedens(step, dgamma2, mu, exp(log_sigma)), id)
  .nll
}
dat <- list(step = trex$step[1:500], id = rep(1:5, each = 100))
par_re <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)),
               log_su = log(0.3), u = rep(0, 5))
obj_re <- MakeADFun(nll_re, par_re, random = "u", silent = TRUE)
sim_re <- simulate(lama(obj_re), random = TRUE, seed = 1)
```
