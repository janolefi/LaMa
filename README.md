
# LaMa - Latent Markov model toolbox 🛠️ <img src="man/figures/Logo_LaMa_surf.png" align="right" style="height:170px; width:auto;">

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/LaMa)](https://CRAN.R-project.org/package=LaMa)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/LaMa)](https://cran.r-project.org/package=LaMa)
[![total
downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/LaMa)](https://cranlogs.r-pkg.org:443/badges/grand-total/LaMa)
[![R-CMD-check](https://github.com/janolefi/LaMa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janolefi/LaMa/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A variety of
<a href="https://journals.sagepub.com/doi/abs/10.1177/1471082X251355681" target="_blank">**latent
Markov models**</a>, including **hidden Markov models** (HMMs), **hidden
semi-Markov models** (HSMMs), **state-space models** (SSMs) and
**continuous-time** variants can be formulated and estimated within the
same framework via directly maximising the likelihood function using the
so-called **forward algorithm**. Applied researchers often need custom
models that standard software does not easily support. Writing tailored
`R` code offers flexibility but suffers from slow estimation speeds.
This `R` package solves these issues by providing easy-to-use
building-blocks for common tasks like the forward algorithm. These
functions can be combined into custom models in a Lego-type approach,
offering very fast estimation speeds. In its most recent iteration,
`LaMa` allows for automatic differentiation with the `RTMB` package
which drastically increases speed and accuracy even more.

<!-- The most important functions are -->

<!-- * `forward` to calculate the log-likelihood for various different models, -->

<!-- * `tpm` for calculating transition probability matrices, -->

<!-- * `stationary` to compute stationary and periodically stationary distributions, -->

<!-- * `stateprobs` and `viterbi` for local and global decoding. -->

## Installation

You can install the released package version from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("LaMa")
```

or the development version from Github:

``` r
remotes::install_github("janolefi/LaMa")
```

<!-- (To install from Github, you need a functional <a href="https://teuder.github.io/rcpp4everyone_en/020_install.html" target="_blank">C++ compiler</a>.) -->

## Package documentation

To aid in building fully custom likelihood functions, this package
contains several vignettes that demonstrate how to simulate data from
and estimate a wide range of models using the functions included in this
package.

HMMs, from simple to complex:

1)  [Introduction to
    LaMa](https://janolefi.github.io/LaMa/articles/Intro_to_LaMa.html)
2)  [Automatic
    differentiation](https://janolefi.github.io/LaMa/articles/LaMa_and_RTMB.html)
3)  [Model
    extensions](https://janolefi.github.io/LaMa/articles/Extensions.html)
4)  [Longitudinal
    data](https://janolefi.github.io/LaMa/articles/Longitudinal_data.html)
5)  [Penalised
    splines](https://janolefi.github.io/LaMa/articles/Penalised_splines.html)

Other latent Markov model classes:

6)  [Continuous-time
    HMMs](https://janolefi.github.io/LaMa/articles/Continuous_time_HMMs.html)
7)  [Markov-modulated Poisson
    processes](https://janolefi.github.io/LaMa/articles/MMMPPs.html)
8)  [State-space
    models](https://janolefi.github.io/LaMa/articles/State_space_models.html)
9)  [Hidden semi-Markov
    models](https://janolefi.github.io/LaMa/articles/HSMMs.html)

<!-- ## Citation -->

<!-- When using LaMa, please cite the package as follows: -->

<!-- ```{r citation} -->

<!-- citation(package = "LaMa") -->

<!-- ``` -->

## Introductory example: Homogeneous HMM

We analyse the `trex` data set contained in the package. It contains
hourly step lengths of a Tyrannosaurus rex, living 66 million years ago.
To these data, we fit a simple 2-state HMM with state-dependent gamma
distributions for the step lengths.

``` r
library(LaMa)
#> Loading required package: RTMB

head(trex, 3)
#>   tod      step     angle state
#> 1   9 0.3252437        NA     1
#> 2  10 0.2458265  2.234562     1
#> 3  11 0.2173252 -2.262418     1
```

We start by defining the negative log-likelihood function. This is made
really convenient by the functions `tpm()` which computes the transition
probability matrix via the multinomial logit link, `stationary()` which
computes the stationary distribution of the Markov chain and `forward()`
which calculates the log-likelihood via the forward algorithm.

``` r
nll = function(par, step){
  # parameter transformations for unconstrained optimisation
  Gamma = tpm(par[1:2]) # rowwise softmax
  delta = stationary(Gamma) # stationary distribution
  mu = exp(par[3:4]) # state-dependent means
  sigma = exp(par[5:6]) # state-dependent sds
  # calculating all state-dependent probabilities
  allprobs = matrix(1, length(step), 2)
  ind = which(!is.na(step))
  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
  # simple forward algorithm to calculate log-likelihood
  -forward(delta, Gamma, allprobs)
}
```

To fit the model, we define the intial parameter vector and numerically
optimise the above function using `nlm()`:

``` r
par = c(-2,-2,             # initial tpm params (logit-scale)
        log(c(0.3, 2.5)),  # initial means for step length (log-transformed)
        log(c(0.2, 1.5)))  # initial sds for step length (log-transformed)

system.time(
  mod <- nlm(nll, par, step = trex$step)
)
#>    user  system elapsed 
#>   0.361   0.012   0.376
```

Really fast for 10.000 data points!

After tranforming the working (unconstrained) parameters to natural
parameters using `tpm()` and `stationary()`, we can visualise the
results:

``` r
# transform parameters to working
(Gamma = tpm(mod$estimate[1:2]))
#>           S1        S2
#> S1 0.8269546 0.1730454
#> S2 0.1608470 0.8391530
(delta = stationary(Gamma)) # stationary HMM
#>       S1       S2 
#> 0.481733 0.518267
(mu = exp(mod$estimate[3:4]))
#> [1] 0.3034926 2.5057054
(sigma = exp(mod$estimate[5:6]))
#> [1] 0.2015258 1.4908153

hist(trex$step, prob = TRUE, bor = "white", breaks = 40, main = "", xlab = "Step length")
curve(delta[1] * dgamma2(x, mu[1], sigma[1]), add = TRUE, lwd = 2, col = "orange", n=500)
curve(delta[2] * dgamma2(x, mu[2], sigma[2]), add = TRUE, lwd = 2, col = "deepskyblue", n=500)
legend("topright", col = c("orange", "deepskyblue"), lwd = 2, bty = "n", legend = c("state 1", "state 2"))
```

<img src="man/figures/README-visualization-1.png" width="75%" style="display: block; margin: auto;" />
