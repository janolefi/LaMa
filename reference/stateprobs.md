# Calculate conditional local state probabilities for homogeneous HMMs

Computes \$\$\Pr(S_t = j \mid X_1, ..., X_T)\$\$ for homogeneous HMMs

## Usage

``` r
stateprobs(
  delta,
  Gamma,
  allprobs,
  trackID = NULL,
  mod = NULL,
  forecast = FALSE
)
```

## Arguments

- delta:

  initial or stationary distribution of length N, or matrix of dimension
  c(k,N) for k independent tracks, if `trackID` is provided

- Gamma:

  transition probability matrix of dimension c(N,N), or array of k
  transition probability matrices of dimension c(N,N,k), if `trackID` is
  provided

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  c(n, N)

- trackID:

  optional vector of length n containing IDs

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. In this case, `Gamma` can be a matrix,
  leading to the same transition probabilities for each track, or an
  array of dimension c(N,N,k), with one (homogeneous) transition
  probability matrix for each track. Furthermore, instead of a single
  vector `delta` corresponding to the initial distribution, a `delta`
  matrix of initial distributions, of dimension c(k,N), can be provided,
  such that each track starts with it's own initial distribution.

- mod:

  optional model object containing initial distribution `delta`,
  transition probability matrix `Gamma`, matrix of state-dependent
  probabilities `allprobs`, and potentially a `trackID` variable

  If you are using automatic differentiation either with
  `RTMB::MakeADFun` or
  [`qreml`](https://janoleko.github.io/reference/qreml.md) and include
  [`forward`](https://janoleko.github.io/reference/forward.md) in your
  likelihood function, the objects needed for state decoding are
  automatically reported after model fitting. Hence, you can pass the
  model object obtained from running `RTMB::report()` or from
  [`qreml`](https://janoleko.github.io/reference/qreml.md) directly to
  this function.

- forecast:

  logical, indicating if forecast probabilities \\\Pr(S_t = j \mid X_1,
  ..., X_t)\\ should be calculated instead.

## Value

matrix of conditional state probabilities of dimension c(n,N)

## See also

Other decoding functions:
[`stateprobs_g()`](https://janoleko.github.io/reference/stateprobs_g.md),
[`stateprobs_p()`](https://janoleko.github.io/reference/stateprobs_p.md),
[`viterbi()`](https://janoleko.github.io/reference/viterbi.md),
[`viterbi_g()`](https://janoleko.github.io/reference/viterbi_g.md),
[`viterbi_p()`](https://janoleko.github.io/reference/viterbi_p.md)

## Examples

``` r
Gamma = tpm(c(-1,-2))
delta = stationary(Gamma)
allprobs = matrix(runif(10), nrow = 10, ncol = 2)

probs = stateprobs(delta, Gamma, allprobs)
```
