# Summary of a `LaMaModel` object

Summarises a model from
[`report`](https://janolefi.github.io/LaMa/reference/report.md) or
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md). Besides
the overview of [`print()`](https://rdrr.io/r/base/print.html), this
contains

- for a homogeneous state process (or one per track) with at most 10
  states: the transition probability matrix (or generator matrix `Q` in
  continuous time, if reported, e.g. by
  [`tpm_ct`](https://janolefi.github.io/LaMa/reference/tpm_ct.md)), its
  stationary distribution and mean dwell times, averaged over tracks if
  there is one per track,

- for a periodic state process: the stationary distribution averaged
  over the cycle, i.e. the long-run proportion of time spent in each
  state,

- the quantities reported by the likelihood function via
  [`REPORT()`](https://rdrr.io/pkg/RTMB/man/TMB-interface.html): numeric
  vectors of length equal to the number of states are shown as one table
  with a column per state, other small objects separately, and large
  ones only by name,

- log-likelihood, AIC, BIC and the maximum gradient component at the
  estimate.

## Usage

``` r
# S3 method for class 'LaMaModel'
summary(object, ...)

# S3 method for class 'summary.LaMaModel'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  model object of class `"LaMaModel"`

- ...:

  additional arguments (not used)

- x:

  object of class `"summary.LaMaModel"`

- digits:

  number of significant digits to print

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.LaMaModel"`, a list containing the quantities shown,
e.g. `Gamma`, `stationary`, `dwell` and `bystate`.
[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
