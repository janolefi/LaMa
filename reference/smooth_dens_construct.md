# Build the design and penalty matrices for smooth density estimation

This high-level function can be used to prepare objects needed to
estimate mixture models of smooth densities using P-Splines.

## Usage

``` r
smooth_dens_construct(
  data,
  par,
  type = "real",
  k = 25,
  knots = NULL,
  degree = 3,
  diff_order = 2
)
```

## Arguments

- data:

  named data frame of 1 or multiple data streams

- par:

  nested named list of initial means and sds/concentrations for each
  data stream

- type:

  vector of length 1 or number of data streams containing the type of
  each data stream, either `"real"` for data on the reals, `"positive"`
  for data on the positive reals or `"circular"` for angular data.

- k:

  vector of length 1 or number of data streams containing the number of
  basis functions for each data stream

- knots:

  optional list of knots vectors (including the boundary knots) to be
  used for basis construction. If not provided, the knots are placed
  equidistantly for `"real"` and `"circular"` and using polynomial
  spacing for `"positive"`.

  For `"real"` and `"positive"` `k - degree + 1` knots are needed, for
  `"circular"` `k + 1` knots are needed.

- degree:

  degree of the B-spline basis functions for each data stream, defaults
  to cubic B-splines

- diff_order:

  order of differencing used for the P-Spline penalty matrix for each
  data stream. Defaults to second-order differences.

## Value

a nested list containing the design matrices `Z`, the penalty matrices
`S`, the initial coefficients `coef` the prediction design matrices
`Z_predict`, the prediction grids `xseq`, and details for the basis
expansion for each data stream.

## Details

Under the hood,
[`make_matrices_dens`](https://janolefi.github.io/LaMa/reference/make_matrices_dens.md)
is used for the actual construction of the design and penalty matrices.

You can provide one or multiple data streams of different types (real,
positive, circular) and specify initial means and standard deviations/
concentrations for each data stream. This information is then converted
into suitable spline coefficients. `smooth_dens_construct` then
constructs the design and penalty matrices for standardised B-splines
basis functions (integrating to one) for each data stream. For types
`"real"` and `"circular"` the knots are placed equidistant in the range
of the data, for type `"positive"` the knots are placed using polynomial
spacing.

## Examples

``` r
## 3 data streams, each with one distribution
# normal data with mean 0 and sd 1
x1 = rnorm(100, mean = 0, sd = 1)
# gamma data with mean 5 and sd 3
x2 = rgamma2(100, mean = 5, sd = 3)
# circular data
x3 = rvm(100, mu = 0, kappa = 2)

data = data.frame(x1 = x1, x2 = x2, x3 = x3)

par = list(x1 = list(mean = 0, sd = 1),
           x2 = list(mean = 5, sd = 3),
           x3 = list(mean = 0, concentration = 2))

SmoothDens = smooth_dens_construct(data, 
                                   par,
                                   type = c("real", "positive", "circular"))
#> x1 
#> Leaving out last column of the penalty matrix, fix the last spline coefficient at zero for identifiability!
#> Parameter matrix excludes the last column. Add a (fixed) zero column using 'cbind(coef, 0)' in your loss function!
#> x2 
#> Leaving out last column of the penalty matrix, fix the last spline coefficient at zero for identifiability!
#> Parameter matrix excludes the last column. Add a (fixed) zero column using 'cbind(coef, 0)' in your loss function!
#> x3 
#> Leaving out last column of the penalty matrix, fix the last spline coefficient at zero for identifiability!
#> Parameter matrix excludes the last column. Add a (fixed) zero column using 'cbind(coef, 0)' in your loss function!
                             
# extracting objects for x1
Z1 = SmoothDens$Z$x1
S1 = SmoothDens$S$x1
coefs1 = SmoothDens$coef$x1

## one data stream, but mixture of two distributions
# normal data with mean 0 and sd 1
x = rnorm(100, mean = 0, sd = 1)
data = data.frame(x = x)

# now parameters for mixture of two normals
par = list(x = list(mean = c(0, 5), sd = c(1,1)))

SmoothDens = smooth_dens_construct(data, par = par)
#> x 
#> Leaving out last column of the penalty matrix, fix the last spline coefficient at zero for identifiability!
#> Parameter matrix excludes the last column. Add a (fixed) zero column using 'cbind(coef, 0)' in your loss function!

# extracting objects 
Z = SmoothDens$Z$x
S = SmoothDens$S$x
coefs = SmoothDens$coef$x
```
