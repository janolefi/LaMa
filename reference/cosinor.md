# Evaluate trigonometric basis expansion

This function can be used to evaluate a trigonometric basis expansion
for a given periodic variable and period. It can also be used in
formulas passed to
[`make_matrices`](https://janolefi.github.io/LaMa/reference/make_matrices.md).

## Usage

``` r
cosinor(x = 1:24, period = 24, eval = TRUE)
```

## Arguments

- x:

  vector of periodic variable values

- period:

  vector of period length. For example for time of day `period = 24`, or
  `period = c(24,12)` for more flexibility.

- eval:

  logical, should not be changed. If `TRUE` the function returns the
  evaluated cosinor terms, if `FALSE` the function returns the terms as
  strings which is used internally form formula evaluation.

## Value

either a desing matrix with the evaluated cosinor terms (`eval = TRUE`)
or a character vector with the terms as strings (`eval = FALSE`).

## Details

The returned basis can be used for linear predictors of the form \$\$
\eta^{(t)} = \beta_0 + \sum\_{k} \bigl( \beta\_{1k} \sin(\frac{2 \pi
t}{period_k}) + \beta\_{2k} \cos(\frac{2 \pi t}{period_k}) \bigr). \$\$
This is relevant for modeling e.g. diurnal variation and the flexibility
can be increased by adding smaller frequencies (i.e. increasing the
length of `period`).

## Examples

``` r
## Evaluate cosinor terms
# builds design matrix
X = cosinor(1:24, period = 24)
X = cosinor(1:24, period = c(24, 12, 6))

## Usage in model formulas
# e.g. frequencies of 24 and 12 hours + interaction with temperature
form = ~ x + temp * cosinor(hour, c(24, 12)) 
data = data.frame(x = runif(24), temp = rnorm(24,20), hour = 1:24)
modmat = make_matrices(form, data = data)
```
