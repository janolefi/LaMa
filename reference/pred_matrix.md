# Build the prediction design matrix based on new data and model_matrices object created by [`make_matrices`](https://janolefi.github.io/reference/make_matrices.md)

Build the prediction design matrix based on new data and model_matrices
object created by
[`make_matrices`](https://janolefi.github.io/reference/make_matrices.md)

## Usage

``` r
pred_matrix(model_matrices, newdata, what = NULL, exclude = NULL)
```

## Arguments

- model_matrices:

  model_matrices object as returned from
  [`make_matrices`](https://janolefi.github.io/reference/make_matrices.md)

- newdata:

  data frame containing the variables in the formula and new data for
  which to evaluate the basis

- what:

  optional character string specifying which formula to use for
  prediction, if `object` contains multiple formulas. If `NULL`, the
  first formula is used.

- exclude:

  optional vector of terms to set to zero in the predicted design
  matrix. Useful for predicting main effects only when e.g.
  `sd(..., bs = "re")` terms are present. See
  [`mgcv::predict.gam`](https://rdrr.io/pkg/mgcv/man/predict.gam.html)
  for more details.

## Value

prediction design matrix for `newdata` with the same basis as used for
`model_matrices`

## Examples

``` r
# single formula
modmat = make_matrices(~ s(x), data.frame(x = 1:10))
Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5))
# with multiple formulas
modmat = make_matrices(list(mu ~ s(x), sigma ~ s(x, bs = "ps")), data = data.frame(x = 1:10))
Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5), what = "mu")
# nested formula list
form = list(stream1 = list(mu ~ s(x), sigma ~ s(x, bs = "ps")))
modmat = make_matrices(form, data = data.frame(x = 1:10))
Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5), what = c("stream1", "mu"))
```
