# Print a `LaMaModel` object

Prints a short overview of a model from
[`lama`](https://janolefi.github.io/LaMa/reference/lama.md) or
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md): the state
process, observations, parameters, log-likelihood, AIC and BIC, and the
names of the reported quantities. Use
[`summary()`](https://rdrr.io/r/base/summary.html) for the estimates.

## Usage

``` r
# S3 method for class 'LaMaModel'
print(x, ...)
```

## Arguments

- x:

  model object of class `"LaMaModel"`

- ...:

  additional arguments (not used)

## Value

`x`, invisibly
