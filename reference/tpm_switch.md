# Build the conditional transition probability matrix of an HSMM

Hidden semi-Markov models are defined in terms of state dwell times and
conditional transition probabilities, given that the current state is
left. These form the matrix \\\Omega\\ (`omega`), also called the
embedded transition probability matrix, whose diagonal entries are zero
as self-transitions are impossible.

This function builds \\\Omega\\ from unconstrained parameters, applying
the inverse multinomial logistic link (softmax) to each row. With a
design matrix `Z`, it builds one matrix for each row of `Z`, so that
\\\Omega\\ depends on covariates. Compatible with automatic
differentiation by `RTMB`.

## Usage

``` r
tpm_switch(beta = NULL, Z = NULL, byrow = FALSE, report = TRUE)
```

## Arguments

- beta:

  parameters; either

  - `NULL` (default) for 2 states, which have no free parameters,

  - a vector of length `nStates * (nStates - 2)`, or

  - a matrix of dimension `c(nStates * (nStates - 2), p + 1)` if `Z` is
    provided, with the intercepts in the first column.

- Z:

  optional covariate design matrix with or without intercept column,
  i.e. of dimension `c(nObs, p)` or `c(nObs, p + 1)`. If it has `p`
  columns, a column of ones is added.

- byrow:

  logical; if `TRUE`, the parameters of each row of the matrix are
  consecutive, see ‘Details’. Defaults to `FALSE`.

- report:

  logical; if `TRUE` (default), the coefficient matrix `beta` is
  reported from the fitted model, if `Z` is provided.

## Value

conditional transition probability matrix of dimension
`c(nStates, nStates)`, or array of such matrices of dimension
`c(nStates, nStates, nObs)` if `Z` is provided

## Details

Each row has `nStates - 1` off-diagonal entries, of which the first is
the reference (its predictor is fixed to 0), so there are
`nStates * (nStates - 2)` free parameters: none for 2 states, 3 for 3
states, 8 for 4 states. The parameters fill a matrix with one row per
state and `nStates - 2` columns, whose row `i` holds the parameters of
the second, third, ... off-diagonal entry of row `i`. By default, this
matrix is filled by column, so that the parameters of row `i` are `i`,
`i + nStates`, ... (as in
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md) of
earlier versions). With `byrow = TRUE`, it is filled by row, so that the
parameters of each row are consecutive. With `Z`, the rows of `beta` are
in the same order.

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`generator_g()`](https://janolefi.github.io/LaMa/reference/generator_g.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_ct()`](https://janolefi.github.io/LaMa/reference/tpm_ct.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# 2 states: no free parameters
omega <- tpm_switch()

# 3 states: 3 free parameters
omega <- tpm_switch(c(0.5, -0.2, 0.1))

# 4 states: 8 free parameters
omega <- tpm_switch(rep(0, 8))

# 3 states, depending on a covariate
z <- rnorm(100)
beta <- cbind(c(0.5, -0.2, 0.1), c(0.3, 0, -0.4)) # intercepts and slopes
omega <- tpm_switch(beta, cbind(1, z)) # array with 100 slices

# 4 states, with the parameters of each row one after another
omega <- tpm_switch(c(0.5, -1, 0.2, 0, 0.3, -0.4, 1, -0.1), byrow = TRUE)
```
