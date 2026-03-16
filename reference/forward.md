# [Forward algorithm](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock) with homogeneous transition probability matrix

Calculates the log-likelihood of a sequence of observations under a
homogeneous hidden Markov model using the **forward algorithm**.

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

  initial or stationary distribution of length `N`, or matrix of
  dimension `c(k,N)` for `k` independent tracks, if `trackID` is
  provided

- Gamma:

  transition probability matrix of dimension `c(N,N)`, or array of `k`
  transition probability matrices of dimension `c(N,N,k)`, if `trackID`
  is provided

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  `c(n, N)`

- trackID:

  optional vector of length `n` containing IDs that separate tracks.

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. In this case, `Gamma` can be a matrix,
  leading to the same transition probabilities for each track, or an
  array of dimension `c(N,N,k)`, with one (homogeneous) transition
  probability matrix for each track. Furthermore, instead of a single
  vector `delta` corresponding to the initial distribution, a `delta`
  matrix of initial distributions, of dimension `c(k,N)`, can be
  provided, such that each track starts with it's own initial
  distribution.

- logspace:

  logical, indicating whether the probabilities/ densities in the
  `allprobs` matrix are on log-scale. If so, internal computations are
  also done on log-scale which is numerically more robust when the
  entries are very small. Note that this is only supported when used in
  AD mode with `RTMB`.

- bw:

  optional integer, indicating the bandwidth for a banded approximation
  of the forward algorithm. This is for expert users only, if sparsity
  in the Hessian matrix w.r.t. observations is required.

- report:

  logical, indicating whether `delta`, `Gamma`, `allprobs`, and
  potentially `trackID` should be reported from the fitted model.
  Defaults to `TRUE`, but only works if `ad = TRUE`, as it uses the
  `RTMB` package.

  When there are multiple tracks, for compatibility with downstream
  functions like
  [`viterbi`](https://janolefi.github.io/LaMa/reference/viterbi.md),
  [`stateprobs`](https://janolefi.github.io/LaMa/reference/stateprobs.md)
  or
  [`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md),
  `forward` should only be called **once** with a `trackID` argument.

- ad:

  optional logical, indicating whether automatic differentiation should
  be used. Determined automatically and intended only for debugging
  purposes.

## Value

HMM log-likelihood for given data and parameters

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
