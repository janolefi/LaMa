# Von Mises distribution

Density, distribution function and random generation for the von Mises
distribution.

## Usage

``` r
dvm(x, mu = 0, kappa = 1, log = FALSE)

pvm(q, mu = 0, kappa = 1, from = NULL, tol = 1e-20)

rvm(n, mu = 0, kappa = 1, wrap = TRUE)
```

## Arguments

- x, q:

  vector of angles measured in radians at which to evaluate the density
  function.

- mu:

  mean direction of the distribution measured in radians.

- kappa:

  non-negative numeric value for the concentration parameter of the
  distribution.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

- from:

  value from which the integration for CDF starts. If `NULL`, is set to
  `mu - pi`.

- tol:

  the precision in evaluating the distribution function

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- wrap:

  logical; if `TRUE`, generated angles are wrapped to the interval
  \[-pi, pi\].

## Value

`dvm` gives the density, `pvm` gives the distribution function, and
`rvm` generates random deviates.

## Details

This implementation of `dvm` allows for automatic differentiation with
`RTMB`. `rvm` and `pvm` are simply wrappers of the corresponding
functions from `circular`.

## Examples

``` r
set.seed(1)
x = rvm(10, 0, 1)
d = dvm(x, 0, 1)
p = pvm(x, 0, 1)
```
