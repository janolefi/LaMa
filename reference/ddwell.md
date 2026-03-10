# State dwell-time distributions of periodically inhomogeneous Markov chains

Computes the dwell-time distribution of a periodically inhomogeneous
Markov chain for a given transition probability matrix.

## Usage

``` r
ddwell(x, Gamma, time = NULL, state = NULL)
```

## Arguments

- x:

  vector of (non-negative) dwell times to compute the dwell-time
  distribution for

- Gamma:

  array of `L` unique transition probability matrices of a periodically
  inhomogeneous Markov chain, with dimensions `c(N,N,L)`, where `N` is
  the number of states and `L` is the cycle length

- time:

  integer vector of time points in `1:L` at which to compute the
  dwell-time distribution. If `NULL`, the overall dwell-time
  distribution is computed.

- state:

  integer vector of state indices for which to compute the dwell-time
  distribution. If `NULL`, dwell-time distributions for all states are
  returned in a named list.

## Value

either time-varying dwell-time distribution(s) if `time` is specified,
or overall dwell-time distribution if `time` is `NULL`. If more than one
`state` is specified, a named list over states is returned.

## Details

For Markov chains whose transition probabilities vary only periodically,
which is achieved for example by expressing the transition probability
matrix as a periodic function of the time of day using
[`tpm_p`](https://janolefi.github.io/LaMa/reference/tpm_p.md) or
[`cosinor`](https://janolefi.github.io/LaMa/reference/cosinor.md), the
probability distribution of time spent in a state can be computed
analytically. This function computes said distribution, either for a
specific time point (conditioning on transitioning into the state at
that time point) or for the overall distribution (conditioning on
transitioning into the state at any time point).

## References

Koslik, J. O., Feldmann, C. C., Mews, S., Michels, R., & Langrock, R.
(2023). Inference on the state process of periodically inhomogeneous
hidden Markov models for animal behavior. arXiv preprint
arXiv:2312.14583.

## Examples

``` r
# setting parameters for trigonometric link
beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
Gamma = tpm_p(beta = beta, degree = 1)

# at specific times and for specific state
ddwell(1:20, Gamma, time = 1:4, state = 1)
#>             1          2         3         4         5          6          7
#> t1 0.06255530 0.08043312 0.1021583 0.1232574 0.1360227 0.13328015 0.11424672
#> t2 0.08580039 0.10897526 0.1314824 0.1450994 0.1421739 0.12187035 0.09171239
#> t3 0.11920292 0.14382237 0.1587175 0.1555173 0.1333083 0.10031987 0.06762680
#> t4 0.16328661 0.18019753 0.1765643 0.1513496 0.1138967 0.07677909 0.04780658
#>             8          9          10          11          12          13
#> t1 0.08597529 0.05795695 0.036086957 0.021451897 0.012526000 0.007332955
#> t2 0.06182440 0.03849502 0.022883373 0.013361855 0.007822280 0.004653506
#> t3 0.04210790 0.02503105 0.014615906 0.008556424 0.005090252 0.003109226
#> t4 0.02841863 0.01659395 0.009714411 0.005779143 0.003530014 0.002237972
#>             14          15           16           17           18           19
#> t1 0.004362405 0.002664643 0.0016893407 0.0011259864 0.0008002580 0.0006147365
#> t2 0.002842453 0.001802070 0.0012011230 0.0008536589 0.0006557576 0.0005503568
#> t3 0.001971199 0.001313852 0.0009337773 0.0007173024 0.0006020094 0.0005554132
#> t4 0.001491663 0.001060150 0.0008143787 0.0006834825 0.0006305802 0.0006396266
#>              20
#> t1 0.0005159290
#> t2 0.0005077585
#> t3 0.0005633812
#> t4 0.0007085138
# results in 4x20 matrix

# or overall distribution for all states
ddwell(1:20, Gamma)
#> $`state 1`
#>          1          2          3          4          5          6          7 
#> 0.23876652 0.16830698 0.11852833 0.08374960 0.05952926 0.04274237 0.03127267 
#>          8          9         10         11         12         13         14 
#> 0.02366237 0.01886914 0.01613048 0.01488364 0.01470017 0.01522215 0.01610698 
#>         15         16         17         18         19         20 
#> 0.01699858 0.01754221 0.01744454 0.01655301 0.01490724 0.01272368 
#> 
#> $`state 2`
#>            1            2            3            4            5            6 
#> 0.5304997868 0.1712799650 0.0662901569 0.0324558381 0.0202817977 0.0156571047 
#>            7            8            9           10           11           12 
#> 0.0141716305 0.0142623506 0.0152462653 0.0166678924 0.0180415973 0.0187787681 
#>           13           14           15           16           17           18 
#> 0.0182904359 0.0162522851 0.0128784300 0.0089320207 0.0053514586 0.0027523266 
#>           19           20 
#> 0.0012168037 0.0004667301 
#> 
# results in list of length 2, each element is a vector of length 20
```
