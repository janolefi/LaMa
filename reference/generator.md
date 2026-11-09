# Build the generator matrix of a continuous-time Markov chain

This function builds the **infinitesimal generator matrix** for a
**continuous-time Markov chain** from an unconstrained parameter vector.

## Usage

``` r
generator(
  beta,
  Z = NULL,
  Eta = NULL,
  byrow = FALSE,
  report = TRUE,
  param = NULL
)
```

## Arguments

- beta:

  parameters; either

  - a vector of length `nStates * (nStates-1)`, or

  - a matrix of dimension `c(nStates * (nStates-1), p+1)` if design
    matrix `Z` is also provided.

- Z:

  optional covariate design matrix with or without intercept column,
  i.e. of dimension `c(nObs, p)` or `c(nObs, p+1)`. If provided, `beta`
  needs to be a matrix of dimension `c(nStates * (nStates-1), p+1)`.

- Eta:

  optional pre-calculated matrix of linear predictors of dimension
  `c(nObs, nStates * (nStates-1))`. If provided, `Z` and `beta` will be
  ignored.

- byrow:

  logical indicating if the generator matrix should be filled by row

- report:

  logical, indicating whether the generator matrix `Q` should be
  reported from the fitted model. Defaults to `TRUE`, but only works if
  when automatic differentiation with `RTMB` is used.

- param:

  depricated, please use argument `beta` instead.

## Value

infinitesimal generator matrix of dimension `c(nStates, nStates)` or
array of such matrices of dimension `c(nStates, nStates, nObs)` if `Z`
or `Eta` is provided.

## Details

**Generator matrices.** The off-diagonal entry \\q\_{ij} \geq 0\\ of a
generator \\Q\\ is the rate of switching from state \\i\\ to state
\\j\\, i.e. \\\Pr(S\_{t+h} = j \mid S_t = i) \approx q\_{ij} h\\ for
small \\h\\. The diagonal entries are \\q\_{ii} = -\sum\_{j \neq i}
q\_{ij}\\, so that each row sums to zero. Hence, the time spent in state
\\i\\ is exponentially distributed with rate \\-q\_{ii}\\, and when
leaving state \\i\\, the chain moves to state \\j\\ with probability
\\q\_{ij} / (-q\_{ii})\\. Transition probabilities over intervals of
length \\\Delta t\\ are \\\exp(Q \Delta t)\\, calculated by
[`tpm_ct`](https://janolefi.github.io/LaMa/reference/tpm_ct.md).

**Parameterisation.** The off-diagonal entries are calculated as
\\q\_{ij} = \exp(\beta_k)\\ to ensure positivity, and the diagonal
entries are then set to the negative row sums. By default
(`byrow = FALSE`), the off-diagonal entries are filled column-wise, i.e.
for 3 states \\\beta\\ corresponds to \\(q\_{21}, q\_{31}, q\_{12},
q\_{32}, q\_{13}, q\_{23})\\. With `byrow = TRUE`, they are filled
row-wise, i.e. \\(q\_{12}, q\_{13}, q\_{21}, q\_{23}, q\_{31},
q\_{32})\\. For 2 states, this means \\\beta = (\log q\_{21}, \log
q\_{12})\\ by default and \\\beta = (\log q\_{12}, \log q\_{21})\\ with
`byrow = TRUE`.

**Covariates.** If a design matrix `Z` is provided, each off-diagonal
entry depends on covariates via \\q\_{ij}^{(t)} = \exp(z_t^\top
\beta\_{ij})\\, where \\z_t\\ is row \\t\\ of `Z`, and one generator
matrix is returned for each row of `Z`. Then `beta` needs to be a matrix
of coefficients of dimension `c(nStates * (nStates-1), p+1)`, where the
first column contains the intercepts and rows are ordered as above.
Alternatively, a matrix of linear predictors `Eta` can be provided, with
one column per off-diagonal entry in the same order. When passed to
[`tpm_ct`](https://janolefi.github.io/LaMa/reference/tpm_ct.md), the
generator is assumed constant within each interval, so row \\t\\ of `Z`
should contain the covariates at the start of interval \\t\\.

## See also

Other transition probability matrix functions:
[`generator_g()`](https://janolefi.github.io/LaMa/reference/generator_g.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_ct()`](https://janolefi.github.io/LaMa/reference/tpm_ct.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# 2 states: 2 free off-diagonal elements
generator(rep(-1, 2))
#>            S1         S2
#> S1 -0.3678794  0.3678794
#> S2  0.3678794 -0.3678794
# 3 states: 6 free off-diagonal elements
generator(rep(-2, 6))
#>            S1         S2         S3
#> S1 -0.2706706  0.1353353  0.1353353
#> S2  0.1353353 -0.2706706  0.1353353
#> S3  0.1353353  0.1353353 -0.2706706
```
