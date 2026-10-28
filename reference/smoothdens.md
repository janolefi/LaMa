# Build the design and penalty matrices for smooth density estimation

This high-level function can be used to prepare objects needed to
estimate mixture models of smooth densities using P-Splines.

## Usage

``` r
smoothdens(
  x,
  par,
  type = "real",
  k = 25,
  knots = NULL,
  degree = 3,
  diff_order = 2,
  pow = 0.5,
  npoints = 10000
)
```

## Arguments

- x:

  data vector for one data stream, or named data frame or list of
  several data streams

- par:

  list of initial means and sds (`"real"`, `"positive"`) or means and
  concentrations (`"circular"`, called `kappa` or `concentration`), one
  value per state. For several data streams, a named list of such lists.

- type:

  type of each data stream, either `"real"` for data on the reals,
  `"positive"` for data on the positive reals or `"circular"` for
  angular data.

- k:

  number of basis functions for each data stream

- knots:

  optional vector of knots (including the boundary knots) to be used for
  basis construction, for several data streams a named list of such
  vectors. If not provided, the knots are placed equidistantly for
  `"real"` and `"circular"` and using polynomial spacing for
  `"positive"`.

  For `"real"` and `"positive"` `k - degree + 1` knots are needed, for
  `"circular"` `k + 1` knots are needed.

- degree:

  degree of the B-spline basis functions for each data stream, defaults
  to cubic B-splines

- diff_order:

  order of differencing used for the P-Spline penalty matrix for each
  data stream. Defaults to second-order differences.

- pow:

  power for polynomial knot spacing, only used for `type = "positive"`

- npoints:

  number of points used in the numerical integration for normalizing the
  B-spline basis functions

  `type`, `k`, `degree`, `diff_order`, `pow` and `npoints` can be of
  length 1 or one per data stream (named or in the order of `x`).

## Value

For one data stream, an object of class `"smoothdens"`: a list
containing the design matrix `Z`, the penalty matrix `S`, the initial
coefficients `coef`, the prediction design matrix `Z_predict`, the
prediction grid `xseq`, details for the basis expansion `basis`, and the
function `density(beta)`. For several data streams, a named list of such
objects.

## Details

Under the hood,
[`make_matrices_dens`](https://janolefi.github.io/LaMa/reference/make_matrices_dens.md)
is used for the actual construction of the design and penalty matrices.

You can provide one or multiple data streams of different types (real,
positive, circular) and specify initial means and standard deviations/
concentrations for each data stream. This information is then converted
into suitable spline coefficients. `smoothdens` then constructs the
design and penalty matrices for standardised B-splines basis functions
(integrating to one) for each data stream. For types `"real"` and
`"circular"` the knots are placed equidistant in the range of the data,
for type `"positive"` the knots are placed using polynomial spacing.

The last spline coefficient of each state is fixed at zero for
identifiability, hence `coef` has `k - 1` columns. `density(beta)` takes
care of this: it maps the coefficient matrix to the `allprobs` matrix,
with rows of missing observations equal to one.

## See also

[`predict`](https://janolefi.github.io/LaMa/reference/predict.smoothdens.md)
for the estimated densities

## Examples

``` r
## one data stream, mixture of two normals
x = c(rnorm(100, 0, 1), rnorm(100, 5, 1))
dens = smoothdens(x, par = list(mean = c(0, 5), sd = c(1, 1)))

Z = dens$Z # design matrix
S = dens$S # penalty matrix
beta = dens$coef # initial coefficients
allprobs = dens$density(beta) # use this in the likelihood function

## 3 data streams
data = data.frame(x1 = rnorm(100, mean = 0, sd = 1),
                  x2 = rgamma2(100, mean = 5, sd = 3),
                  x3 = rvm(100, mu = 0, kappa = 2))

par = list(x1 = list(mean = 0, sd = 1),
           x2 = list(mean = 5, sd = 3),
           x3 = list(mean = 0, kappa = 2))

dens = smoothdens(data, par, type = c("real", "positive", "circular"))

# objects for x1
Z1 = dens$x1$Z
allprobs = dens$x1$density(dens$x1$coef) *
  dens$x2$density(dens$x2$coef) *
  dens$x3$density(dens$x3$coef)
```
