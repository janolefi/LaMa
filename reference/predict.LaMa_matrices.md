# Build the prediction design matrix based on new data and model_matrices object created by [`make_matrices`](https://janolefi.github.io/LaMa/reference/make_matrices.md)

Build the prediction design matrix based on new data and model_matrices
object created by
[`make_matrices`](https://janolefi.github.io/LaMa/reference/make_matrices.md)

## Usage

``` r
# S3 method for class 'LaMa_matrices'
predict(object, newdata, what = NULL, ...)
```

## Arguments

- object:

  model matrices object as returned from
  [`make_matrices`](https://janolefi.github.io/LaMa/reference/make_matrices.md)

- newdata:

  data frame containing the variables in the formula and new data for
  which to evaluate the basis

- what:

  optional character string specifying which formula to use for
  prediction if `object` contains multiple formulas.

- ...:

  for method consistency only

## Value

prediction design matrix for `newdata` with the same basis as used for
`model_matrices`

## See also

[`make_matrices`](https://janolefi.github.io/LaMa/reference/make_matrices.md)
for creating objects of class `LaMa_matrices` which can be used for
prediction by this function.

## Examples

``` r
# single formula
modmat = make_matrices(~ s(x), data.frame(x = 1:10))
Z_p = predict(modmat, data.frame(x = 1:10 - 0.5))
# with multiple formulas
modmat = make_matrices(list(mu ~ s(x), sigma ~ s(x, bs = "ps")), data = data.frame(x = 1:10))
Z_p = predict(modmat, data.frame(x = 1:10 - 0.5), what = "mu")
# nested formula list
form = list(stream1 = list(mu ~ s(x), sigma ~ s(x, bs = "ps")))
modmat = make_matrices(form, data = data.frame(x = 1:10))
Z_p = predict(modmat, data.frame(x = 1:10 - 0.5), what = c("stream1", "mu"))
```
