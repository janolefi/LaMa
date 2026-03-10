# Build all transition probability matrices of an inhomogeneous HMM

In an HMM, we often model the influence of covariates on the state
process by linking them to the transition probabiltiy matrix. Most
commonly, this is done by specifying a linear predictor \$\$
\eta\_{ij}^{(t)} = \beta^{(ij)}\_0 + \beta^{(ij)}\_1 z\_{t1} + \dots +
\beta^{(ij)}\_p z\_{tp} \$\$ for each off-diagonal element (\\i \neq
j\\) of the transition probability matrix and then applying the inverse
multinomial logistic link (also known as softmax) to each row. This
function efficiently calculates all transition probabilty matrices for a
given design matrix `Z` and parameter matrix `beta`.

## Usage

``` r
tpm_g(
  Z,
  beta,
  Eta = NULL,
  byrow = FALSE,
  ref = NULL,
  ad = NULL,
  report = TRUE,
  sparse = FALSE
)
```

## Arguments

- Z:

  covariate design matrix with or without intercept column, i.e. of
  dimension c(n, p) or c(n, p+1)

  If `Z` has only p columns, an intercept column of ones will be added
  automatically.

- beta:

  matrix of coefficients for the off-diagonal elements of the transition
  probability matrix

  Needs to be of dimension c(N\*(N-1), p+1), where the first column
  contains the intercepts.

- Eta:

  optional pre-calculated linear predictor matrix of dimension c(n,
  N\*(N-1)).

  Usually, `Eta` is calculated as `Z %*% t(beta)`. If provided, no `Z`
  and `beta` are necessary and will be ignored.

- byrow:

  logical indicating if each transition probability matrix should be
  filled by row

  Defaults to `FALSE`, but should be set to `TRUE` if one wants to work
  with a matrix of beta parameters returned by popular HMM packages like
  `moveHMM`, `momentuHMM`, or `hmmTMB`.

- ref:

  Optional integer vector of length N giving, for each row, the column
  index of the reference state (its predictor is fixed to 0). Defaults
  to the diagonal (`ref = 1:N`).

- ad:

  optional logical, indicating whether automatic differentiation with
  `RTMB` should be used. By default, the function determines this
  itself.

- report:

  logical, indicating whether the coefficient matrix `beta` should be
  reported from the fitted model. Defaults to `TRUE`, but only works if
  `ad = TRUE`.

- sparse:

  logical, indicating whether sparsity in the rows of `Z` should be
  exploited.

## Value

array of transition probability matrices of dimension c(N,N,n)

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_cont()`](https://janolefi.github.io/LaMa/reference/tpm_cont.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
Z = matrix(runif(200), ncol = 2)
beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
Gamma = tpm_g(Z, beta)
```
