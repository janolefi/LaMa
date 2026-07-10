# Skew normal distribution

Density, distribution function, quantile function and random generation
for the skew normal distribution.

## Usage

``` r
dskewnorm(x, xi = 0, omega = 1, alpha = 0, log = FALSE)

pskewnorm(q, xi = 0, omega = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE)

qskewnorm(p, xi = 0, omega = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE)

rskewnorm(n, xi = 0, omega = 1, alpha = 0)
```

## Arguments

- x, q:

  vector of quantiles

- xi:

  location parameter

- omega:

  scale parameter, must be positive.

- alpha:

  skewness parameter, +/- `Inf` is allowed.

- log, log.p:

  logical; if `TRUE`, probabilities/ densities \\p\\ are returned as
  \\\log(p)\\.

- lower.tail:

  logical; if `TRUE` (default), probabilities are \\P\[X \le x\]\\,
  otherwise \\P\[X \> x\]\\.

- p:

  vector of probabilities

- n:

  number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

## Value

`dskewnorm` gives the density, `pskewnorm` gives the distribution
function, `qskewnorm` gives the quantile function, and `rskewnorm`
generates random deviates.

## Details

This implementation of `dskewnorm` allows for automatic differentiation
with `RTMB` while the other functions are imported from the `sn`
package.

## Examples

``` r
x = rskewnorm(1)
d = dskewnorm(x)
p = pskewnorm(x)
q = qskewnorm(p)
```
