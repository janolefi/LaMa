# Sparsity-retaining matrix multiplication

Standard matrix multiplication destroys automatic sparsity detection by
`RTMB` which is essential for models with high-dimensional random
effects. This can be mitigated by changing to "plain" with `TapeConfig`,
but this can make AD tape construction very slow. Here, we provide a
different version that retains sparsity. It may be slightly slower than
the standard method when constructing the AD tape.

## Usage

``` r
A %sp% B
```

## Arguments

- A:

  matrix of dimension `c(n, p)`

- B:

  matrix of dimension `c(p, m)`

## Value

the matrix product of A and B, which is of dimension `c(n, m)`

## Examples

``` r
A <- matrix(1:6, nrow = 2, ncol = 3)
B <- matrix(7:12, nrow = 3, ncol = 2)
A %sp% B
#>      [,1] [,2]
#> [1,]   76  103
#> [2,]  100  136
```
