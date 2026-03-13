# Sample from a multivariate Gaussian with a sparse precision matrix

Draw samples from a multivariate Gaussian distribution specified by a
sparse precision matrix. This is numerically efficient for
high-dimensional but sparse systems.

## Usage

``` r
rgmrf(n, mean = 0, Q)
```

## Arguments

- n:

  Number of samples to draw.

- mean:

  Mean vector (or scalar, which will be recycled to match the dimension
  of `Q`).

- Q:

  Sparse precision matrix (\\\Sigma^{-1}\\).

## Value

A matrix of samples with rows corresponding to samples and columns to
dimensions.

## Examples

``` r
rgmrf(3, mean = c(1, 2, 3), Q = Matrix::Diagonal(3))
#>          [,1]       [,2]     [,3]
#> [1,] 1.760995 -0.3062566 4.064045
#> [2,] 1.343252  1.5453050 3.537838
#> [3,] 1.279675  2.6315274 3.529627
```
