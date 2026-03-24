# Forward algorithm to calculate the HMM log-likelihood

Calculates the log-likelihood of a sequence of observations under a
hidden Markov model using the **forward algorithm** (Zucchini, MacDonald
& Langrock, 2016).

## Usage

``` r
forward(
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

  transition probability matrix; either

  - a matrix of dimension `c(nStates, nStates)`, or

  - an array of dimension `c(nStates, nStates, nTracks)` if `trackID` is
    provided, or

  - an array of dimension `c(nStates, nStates, nObs)` for time-varying
    transition probabilities, in which case
    [`forward_g`](https://janolefi.github.io/LaMa/reference/forward_g.md)
    is called internally.

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

If `trackID` is provided, the total log-likelihood is the sum of each
track's likelihood contribution. In this case, `Gamma` can be

- a matrix (same transition probabilities for each track),

- an array of dimension `c(nStates, nStates, nTracks)` (track-specific
  transition probability matrix), or

- an array of dimension `c(nStates, nStates, nObs)` for time-varying
  transition probabilities, in which case
  [`forward_g`](https://janolefi.github.io/LaMa/reference/forward_g.md)
  is called internally.

Additionally, `delta` can be a vector (same initial distribution for
each track) or a matrix of dimension `c(nTracks, nStates)`
(track-specific initial distributions).

**Note:** When there are multiple tracks, for compatibility with
downstream functions like
[`viterbi`](https://janolefi.github.io/LaMa/reference/viterbi.md),
[`stateprobs`](https://janolefi.github.io/LaMa/reference/stateprobs.md)
or
[`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md),
`forward` should only be called **once** with a `trackID` argument.

## References

Zucchini, W., MacDonald, I.L., & Langrock, R. (2016). *Hidden Markov
Models for Time Series: An Introduction Using R* (2nd ed.). Chapman &
Hall/CRC.

## See also

Other forward algorithms:
[`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md),
[`forward_hsmm()`](https://janolefi.github.io/LaMa/reference/forward_hsmm.md),
[`forward_ihsmm()`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md),
[`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md),
[`forward_phsmm()`](https://janolefi.github.io/LaMa/reference/forward_phsmm.md)

## Examples

``` r
## negative log likelihood function
nll = function(par, step) {
 # parameter transformations for unconstrained optimisation
 Gamma = tpm(par[1:2]) # multinomial logit link
 delta = stationary(Gamma) # stationary HMM
 mu = exp(par[3:4])
 sigma = exp(par[5:6])
 # calculate all state-dependent probabilities
 allprobs = matrix(1, length(step), 2)
 ind = which(!is.na(step))
 for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
 # simple forward algorithm to calculate log-likelihood
 -forward(delta, Gamma, allprobs)
}

## fitting an HMM to the trex data
par = c(-2,-2,            # initial tpm params (logit-scale)
        log(c(0.3, 2.5)), # initial means for step length (log-transformed)
        log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
mod = nlm(nll, par, step = trex$step[1:1000])
```
