# State-dependent densities of one data stream

Builds the matrix of state-dependent densities (or probabilities) of one
data stream, as used in the `allprobs` matrix of
[`forward`](https://janolefi.github.io/LaMa/reference/forward.md).
Several data streams that are conditionally independent given the state
are combined by multiplying their matrices.

## Usage

``` r
statedens(x, dens, ..., nStates = NULL, log = FALSE)
```

## Arguments

- x:

  vector of observations of length `nObs`, possibly containing `NA`s

- dens:

  density function, e.g. `dnorm`,
  [`dgamma2`](https://janolefi.github.io/LaMa/reference/gamma2.md) or
  [`dvm`](https://janolefi.github.io/LaMa/reference/vm.md), with the
  observations as first argument

- ...:

  parameters of `dens`, see ‘Details’

- nStates:

  number of states. Only needed if all parameters are single values.

- log:

  logical; if `TRUE`, log-densities are returned

## Value

matrix of dimension `c(nObs, nStates)`

## Details

Each parameter in `...` can be

- a single value, used for all states,

- a vector of length `nStates`, one value per state, or

- a matrix of dimension `c(nObs, nStates)`, one value per observation
  and state, e.g. from a regression on covariates.

Missing observations contribute nothing to the likelihood, so their rows
are 1 (0 with `log = TRUE`). As each data stream gets its own matrix, a
missing value in one stream does not remove the observations of the
others.

With `log = TRUE`, the log-densities are returned, to be added over data
streams and used with `forward(..., logspace = TRUE)`. They are computed
by `dens(..., log = TRUE)` if `dens` has an argument `log`, otherwise as
`log(dens(...))` (e.g. for densities from
[`zero_inflate`](https://janolefi.github.io/LaMa/reference/zero_inflate.md)).

In `obj$simulate()` of `RTMB`, `x` is registered as data and simulated
by [`forward`](https://janolefi.github.io/LaMa/reference/forward.md),
see
[`simulate.LaMaModel`](https://janolefi.github.io/LaMa/reference/simulate.LaMaModel.md).

Writing the loop over states by hand is equally valid and gives the same
result. `statedens()` is just a shortcut.

## Examples

``` r
step <- trex$step[1:200]
angle <- trex$angle[1:200]
mu <- c(0.3, 2.5); sigma <- c(0.2, 1.5); kappa <- c(0.2, 1.5)

# one data stream
allprobs <- statedens(step, dgamma2, mu, sigma)

# the same by hand
allprobs2 <- matrix(1, length(step), 2)
ind <- which(!is.na(step))
for(j in 1:2) allprobs2[ind, j] <- dgamma2(step[ind], mu[j], sigma[j])
all.equal(allprobs, allprobs2)
#> [1] TRUE

# two data streams, mean of the turning angle the same in both states
allprobs <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)

# on the log scale, for forward(..., logspace = TRUE)
logprobs <- statedens(step, dgamma2, mu, sigma, log = TRUE) +
  statedens(angle, dvm, 0, kappa, log = TRUE)
```
