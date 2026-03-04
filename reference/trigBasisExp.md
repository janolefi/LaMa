# Compute the design matrix for a trigonometric basis expansion

Given a periodically varying variable such as time of day or day of year
and the associated cycle length, this function performs a basis
expansion to efficiently calculate a linear predictor of the form \$\$
\eta^{(t)} = \beta_0 + \sum\_{k=1}^K \bigl( \beta\_{1k} \sin(\frac{2 \pi
k t}{L}) + \beta\_{2k} \cos(\frac{2 \pi k t}{L}) \bigr). \$\$ This is
relevant for modeling e.g. diurnal variation and the flexibility can be
increased by adding smaller frequencies (i.e. increasing \\K\\).

## Usage

``` r
trigBasisExp(tod, L = 24, degree = 1)
```

## Arguments

- tod:

  equidistant sequence of a cyclic variable

  For time of day and e.g. half-hourly data, this could be 1, ..., L and
  L = 48, or 0.5, 1, 1.5, ..., 24 and L = 24.

- L:

  length of one cycle on the scale of the time variable. For time of
  day, this would be 24.

- degree:

  degree K of the trigonometric link above. Increasing K increases the
  flexibility.

## Value

design matrix (without intercept column), ordered as sin1, cos1, sin2,
cos2, ...

## Examples

``` r
# no examples
```
