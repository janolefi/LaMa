# [Forward algorithm](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock) for homogeneous hidden semi-Markov models

Calculates the (approximate) log-likelihood of a sequence of
observations under a homogeneous hidden semi-Markov model using a
modified **forward algorithm**.

## Usage

``` r
forward_hsmm(
  dm,
  omega,
  allprobs,
  trackID = NULL,
  delta = NULL,
  eps = 1e-10,
  report = TRUE
)
```

## Arguments

- dm:

  list of length N containing vectors of dwell-time probability mass
  functions (PMFs) for each state. The vector lengths correspond to the
  approximating state aggregate sizes, hence there should be little
  probablity mass not covered by these.

- omega:

  matrix of dimension c(N,N) of conditional transition probabilites,
  also called embedded transition probability matrix.

  Contains the transition probabilities given that the current state is
  left. Hence, the diagonal elements need to be zero and the rows need
  to sum to one. Can be constructed using
  [`tpm_emb`](https://janolefi.github.io/LaMa/reference/tpm_emb.md).

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  c(n, N) which will automatically be converted to the appropriate
  dimension.

- trackID:

  optional vector of length n containing IDs

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. In this case, `dm` can be a nested list,
  where the top layer contains k `dm` lists as described above. `omega`
  can then also be an array of dimension c(N,N,k) with one conditional
  transition probability matrix for each track. Furthermore, instead of
  a single vector `delta` corresponding to the initial distribution, a
  `delta` matrix of initial distributions, of dimension c(k,N), can be
  provided, such that each track starts with it's own initial
  distribution.

- delta:

  optional vector of initial state probabilities of length N

  By default, the stationary distribution is computed (which is
  typically recommended).

- eps:

  small value to avoid numerical issues in the approximating transition
  matrix construction. Usually, this should not be changed.

- report:

  logical, indicating whether initial distribution, approximating
  transition probability matrix and `allprobs` matrix should be reported
  from the fitted model. Defaults to `TRUE`.

## Value

HSMM log-likelihood for given data and parameters

## Details

Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs,
where the state duration distribution is explicitly modelled by a
distribution on the positive integers. For direct numerical maximum
likelhood estimation, HSMMs can be represented as HMMs on an enlarged
state space (of size \\M\\) and with structured transition
probabilities.

This function is designed to be used with automatic differentiation
based on the `R` package `RTMB`. It will be very slow without it!

## References

Langrock, R., & Zucchini, W. (2011). Hidden Markov models with arbitrary
state dwell-time distributions. Computational Statistics & Data
Analysis, 55(1), 715-724.

Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state
dwell-time distributions. Computational Statistics & Data Analysis, 209,
108171.

## See also

Other forward algorithms:
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
[`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md),
[`forward_ihsmm()`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md),
[`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md),
[`forward_phsmm()`](https://janolefi.github.io/LaMa/reference/forward_phsmm.md)

## Examples

``` r
# currently no examples
```
