# Build the generator matrix of a continuous-time Markov chain

This function builds the **infinitesimal generator matrix** for a
**continuous-time Markov chain** from an unconstrained parameter vector.

## Usage

``` r
generator(param, byrow = FALSE, report = TRUE)
```

## Arguments

- param:

  unconstrained parameter vector of length N\*(N-1) where N is the
  number of states of the Markov chain

- byrow:

  logical indicating if the transition probability matrix should be
  filled by row

- report:

  logical, indicating whether the generator matrix Q should be reported
  from the fitted model. Defaults to `TRUE`, but only works if when
  automatic differentiation with `RTMB` is used.

## Value

infinitesimal generator matrix of dimension c(N,N)

## See also

Other transition probability matrix functions:
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_cont()`](https://janolefi.github.io/LaMa/reference/tpm_cont.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# 2 states: 2 free off-diagonal elements
generator(rep(-1, 2))
#>            S1         S2
#> S1 -0.3678794  0.3678794
#> S2  0.3678794 -0.3678794
# 3 states: 6 free off-diagonal elements
generator(rep(-2, 6))
#>            S1         S2         S3
#> S1 -0.2706706  0.1353353  0.1353353
#> S2  0.1353353 -0.2706706  0.1353353
#> S3  0.1353353  0.1353353 -0.2706706
```
