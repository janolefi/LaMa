# Calculate continuous time transition probabilities

A continuous-time Markov chain is described by an infinitesimal
generator matrix \\Q\\. When observing data at time points \\t_1, \dots,
t_n\\ the transition probabilities between \\t_i\\ and \\t\_{i+1}\\ are
calculated as \$\$\Gamma(\Delta t_i) = \exp(Q \Delta t_i),\$\$ where
\\\exp()\\ is the matrix exponential. The mapping \\\Gamma(\Delta t)\\
is also called the **Markov semigroup**. This function calculates all
transition matrices based on a given generator and time differences.

## Usage

``` r
tpm_ct(Q, timediff, rates = NULL, ad = NULL, report = TRUE)

tpm_cont(Q, timediff, rates = NULL, ad = NULL, report = TRUE)
```

## Arguments

- Q:

  infinitesimal generator matrix of the continuous-time Markov chain of
  dimension `c(nStates, nStates)`, or array of generator matrices of
  dimension `c(nStates, nStates, nObs-1)` (or larger) for time-varying
  generators, e.g. from
  [`generator`](https://janolefi.github.io/LaMa/reference/generator.md)`(beta, Z)`.

- timediff:

  time differences between observations of length `nObs`-1 when based on
  `nObs` observations

- rates:

  optional vector of state-dependent rates for MM(M)PP fitting. For the
  MM(M)PP likelihood, the matrices needed in the forward algorithm are
  \\\exp((Q - \Lambda) \Delta t)\\, where \\\Lambda\\ is a diagonal
  matrix with the state-dependent rates on the diagonal.

- ad:

  optional logical, indicating whether automatic differentiation with
  `RTMB` should be used. By default, the function determines this
  itself.

- report:

  logical, indicating whether `Q` should be reported from the fitted
  model. Defaults to `TRUE`, but only works if `ad = TRUE`.

## Value

array of continuous-time transition matrices of dimension
`c(nStates, nStates, nObs-1)`

## Details

**Transition probabilities.** The entry \\(i,j)\\ of \\\Gamma(\Delta t)
= \exp(Q \Delta t)\\ is the probability that the chain is in state \\j\\
at time \\t + \Delta t\\, given that it is in state \\i\\ at time \\t\\.
Its rows sum to one. For short intervals, \\\Gamma(\Delta t)\\ is close
to the identity matrix. For long intervals, each row is close to the
stationary distribution (see
[`stationary_ct`](https://janolefi.github.io/LaMa/reference/stationary_ct.md)).
This is what makes continuous-time models suitable for irregularly
sampled data: each interval gets its own transition matrix.

**Time-varying generators.** If `Q` is an array of generator matrices,
e.g. from
[`generator`](https://janolefi.github.io/LaMa/reference/generator.md)`(beta, Z)`,
the generator is assumed constant within each interval, and `Q[, , t]`
is used for `timediff[t]`. Hence, row \\t\\ of `Z` should contain the
covariates at the start of interval \\t\\, i.e. at observation \\t\\.
Only the first `length(timediff)` slices of `Q` are used.

**Repeated time differences.** If `Q` is a single matrix,
\\\Gamma(\Delta t)\\ only depends on \\\Delta t\\. Hence, the matrix
exponential is computed only once for each unique value in `timediff`
and then copied to all intervals with that time difference. This happens
automatically, so there is no need to do it manually. It saves a lot of
time when there are only few unique time differences, e.g. when time is
measured in days or observations are regular with occasional gaps.
Copying is free on the `RTMB` tape, and when all time differences are
unique, nothing changes.

**MM(M)PPs.** If `rates` are provided, the function returns \\\exp((Q -
\Lambda) \Delta t)\\ with \\\Lambda = \mathrm{diag}(\lambda_1, \dots,
\lambda_N)\\. Its entry \\(i,j)\\ is the probability of moving from
state \\i\\ to state \\j\\ over the interval without any arrival in
between, so its rows sum to less than one. Multiplied with the
state-dependent rates at the arrival times (in `allprobs`), these
matrices give the MM(M)PP likelihood in the forward algorithm.

**Computation.** Without automatic differentiation, the matrix
exponential is calculated in C++ with Armadillo, separately for each
interval. With automatic differentiation (`ad = TRUE`),
[`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html) is
called for each interval, i.e. a Padé approximation with scaling and
squaring, all of which is recorded on the `RTMB` tape.

For **two states**, a closed form is used in the AD case instead. Let
\\A = Q - \Lambda\\ (or \\A = Q\\ without rates) with real eigenvalues
\\\nu_1 \> \nu_2\\. Then \$\$\exp(A \Delta t) = e^{\nu_1 \Delta t} I +
\frac{e^{\nu_1 \Delta t} - e^{\nu_2 \Delta t}}{\nu_1 - \nu_2} (A - \nu_1
I).\$\$ Without rates, \\\nu_1 = 0\\ and \\\nu_2 = -(q\_{12} +
q\_{21})\\, so that, e.g., \\\Gamma\_{12}(\Delta t) = q\_{12} (1 -
e^{-(q\_{12} + q\_{21}) \Delta t}) / (q\_{12} + q\_{21})\\. This is
exact, not an approximation. For two states, models are therefore
expected to be much faster to fit (often by more than an order of
magnitude for the gradient): the closed form needs only a few
elementwise operations on vectors of length `length(timediff)`, whereas
[`Matrix::expm`](https://rdrr.io/pkg/Matrix/man/expm-methods.html)
records a full Padé approximation (matrix products, a linear solve and
repeated squaring) on the tape separately for each interval, all of
which needs to be evaluated and differentiated in every optimisation
step.

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`generator_g()`](https://janolefi.github.io/LaMa/reference/generator_g.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# building a Q matrix for a 3-state cont.-time Markov chain
Q = generator(rep(-2, 6))

# draw random time differences
timediff = rexp(100, 10)

# compute all transition matrices
Gamma = tpm_cont(Q, timediff)
#> tpm_cont() is deprecated, please use tpm_ct() instead.
```
