# Forward algorithm with time-varying transition probability matrix

Calculates the log-likelihood of a sequence of observations under a
hidden Markov model with time-varying transition probabilities using the
**forward algorithm** (Zucchini, MacDonald & Langrock, 2016).

## Usage

``` r
forward_g(
  delta,
  Gamma,
  allprobs,
  trackID = NULL,
  logspace = FALSE,
  bw = NULL,
  report = TRUE,
  ad = NULL
)
```

## Arguments

- delta:

  initial distribution; either

  - a vector of length `nStates`, or

  - a matrix of dimension `c(nTracks, nStates)`, if `trackID` is
    provided.

- Gamma:

  array of transition probability matrices of dimension
  `c(nStates, nStates, nObs)`, where the first slice of each track is
  ignored as there is no transition into the start of a track.

  For a single track, an array of dimension
  `c(nStates, nStates, nObs-1)` is also accepted.

  This function also supports continuous-time HMMs, where each slice is
  a Markov semigroup \\\Gamma(\Delta t) = \exp(Q \Delta t)\\ for
  generator \\Q\\.

- allprobs:

  matrix of state-dependent probabilities or density values of dimension
  `c(nObs, nStates)`

- trackID:

  optional vector of length `nObs` containing `nTracks` unique IDs that
  separate tracks (see ‘Details’).

- logspace:

  logical; if `TRUE`, `allprobs` is assumed to be on the log-scale,
  improving numerical stability for small probabilities. Only supported
  with `RTMB`.

- bw:

  optional positive integer specifying the bandwidth for a banded
  approximation of the forward algorithm, inducing a banded Hessian
  w.r.t. the observations. Defaults to `NULL` (exact algorithm).
  Approximation error decays geometrically in `bw`.

- report:

  logical; if `TRUE` (default), `delta`, `Gamma`, `allprobs`, and
  `trackID` are reported from the fitted model. Requires `ad = TRUE`.

- ad:

  logical; whether to use automatic differentiation. Determined
  automatically — for debugging only.

## Value

HMM log-likelihood for given data and parameters

## Details

If `trackID` is provided, the total log-likelihood will be the sum of
each track's likelihood contribution. In this case, `Gamma` must be an
array of dimension `c(nStates, nStates, nObs)`, matching the number of
rows of allprobs. For each track, the transition matrix at the beginning
of the track will be ignored (as there is no transition between tracks).
Additionally, `delta` can be a vector (same initial distribution for
each track) or a matrix of dimension `c(nTracks, nStates)` (different
initial distribution for each track).

**Note:** When there are multiple tracks, for compatibility with
downstream functions like
[`viterbi_g`](https://janolefi.github.io/LaMa/reference/viterbi_g.md),
[`stateprobs_g`](https://janolefi.github.io/LaMa/reference/stateprobs_g.md)
or
[`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md),
`forward_g` should only be called **once** with a `trackID` argument.

## References

Zucchini, W., MacDonald, I.L., & Langrock, R. (2016). *Hidden Markov
Models for Time Series: An Introduction Using R* (2nd ed.). Chapman &
Hall/CRC.

Mews, Koslik, & Langrock (2025). How to build your latent Markov model:
The role of time and space. Statistical Modelling, 25(6), 481–507.

## See also

Other forward algorithms:
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
[`forward_hsmm()`](https://janolefi.github.io/LaMa/reference/forward_hsmm.md),
[`forward_ihsmm()`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md),
[`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md),
[`forward_phsmm()`](https://janolefi.github.io/LaMa/reference/forward_phsmm.md)

## Examples

``` r
## Simple usage
Gamma = array(c(0.9, 0.2, 0.1, 0.8), dim = c(2,2,10))
delta = c(0.5, 0.5)
allprobs = matrix(0.5, 10, 2)
forward_g(delta, Gamma, allprobs)
#> [1] -6.931472
# \donttest{
## Full model fitting example
## negative log likelihood function
nll = function(par, step, Z) {
 # parameter transformations for unconstrained optimisation
 beta = matrix(par[1:6], nrow = 2)
 Gamma = tpm_g(Z, beta) # multinomial logit link for each time point
 delta = stationary(Gamma[,,1]) # stationary HMM
 mu = exp(par[7:8])
 sigma = exp(par[9:10])
 # calculate all state-dependent probabilities
 allprobs = matrix(1, length(step), 2)
 ind = which(!is.na(step))
 for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
 # simple forward algorithm to calculate log-likelihood
 -forward_g(delta, Gamma, allprobs)
}

## fitting an HMM to the trex data
par = c(-1.5,-1.5,        # initial tpm intercepts (logit-scale)
        rep(0, 4),        # initial tpm slopes
        log(c(0.3, 2.5)), # initial means for step length (log-transformed)
        log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
mod = nlm(nll, par, step = trex$step[1:500], Z = cosinor(trex$tod[1:500]))
#> Warning: NA/NaN replaced by maximum positive value
# }
```
