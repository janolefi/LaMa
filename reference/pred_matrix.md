# Build the prediction design matrix based on new data and model_matrices object created by [`modmat`](https://janolefi.github.io/LaMa/reference/modmat.md)

**Deprecated.** Use
[`predict`](https://janolefi.github.io/LaMa/reference/predict.LaMa_matrices.md)
on the object returned by
[`modmat`](https://janolefi.github.io/LaMa/reference/modmat.md) instead,
which does the same thing: `predict(modmat, newdata)` in place of
`pred_matrix(modmat, newdata)`.

## Usage

``` r
pred_matrix(model_matrices, newdata, what = NULL, exclude = NULL)
```

## Arguments

- model_matrices:

  model_matrices object as returned from
  [`modmat`](https://janolefi.github.io/LaMa/reference/modmat.md)

- newdata:

  data frame containing the variables in the formula and new data for
  which to evaluate the basis

- what:

  optional character string specifying which formula to use for
  prediction if `model_matrices` contains multiple formulas.

- exclude:

  optional vector of terms to set to zero in the predicted design
  matrix. Useful for predicting main effects only when e.g.
  `sd(..., bs = "re")` terms are present. See
  [`mgcv::predict.gam`](https://rdrr.io/pkg/mgcv/man/predict.gam.html)
  for more details.

## Value

prediction design matrix for `newdata` with the same basis as used for
`model_matrices`

## See also

[`predict`](https://janolefi.github.io/LaMa/reference/predict.LaMa_matrices.md),
which replaces this function

## Examples

``` r
modmat = modmat(~ s(x), data.frame(x = 1:10))
# deprecated:
# Z_p = pred_matrix(modmat, data.frame(x = 1:10 - 0.5))
# use instead:
Z_p = predict(modmat, data.frame(x = 1:10 - 0.5))
```
