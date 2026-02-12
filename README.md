
# LaMa - Latent Markov model toolbox 🛠️ <img src="man/figures/Logo_LaMa_surf.png" align="right" style="height:170px; width:auto;">

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/LaMa)](https://CRAN.R-project.org/package=LaMa)
[![metacran
downloads](https://cranlogs.r-pkg.org/badges/last-month/LaMa)](https://cran.r-project.org/package=LaMa)
[![total
downloads](https://cranlogs.r-pkg.org:443/badges/grand-total/LaMa)](https://cranlogs.r-pkg.org:443/badges/grand-total/LaMa)
[![R-CMD-check](https://github.com/janoleko/LaMa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/janoleko/LaMa/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

A variety of **latent Markov models**
<a href="https://journals.sagepub.com/doi/abs/10.1177/1471082X251355681" target="_blank">(Mews,
Koslik, and Langrock 2025)</a>, including **hidden Markov models**
(HMMs), **hidden semi-Markov models** (HSMMs), **state-space models**
(SSMs) and **continuous-time** variants can be formulated and estimated
within the same framework via directly maximising the likelihood
function using the so-called **forward algorithm**
<a href="https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock" target="_blank">(Zucchini,
MacDonald, and Langrock 2016)</a>. Applied researchers often need custom
models that standard software does not easily support. Writing tailored
`R` code offers flexibility but suffers from slow estimation speeds.
This `R` package solves these issues by providing easy-to-use functions
(written in C++ for speed) for common tasks like the forward algorithm.
These functions can be combined into custom models in a Lego-type
approach, offering up to 10-20 times faster estimation via standard
numerical optimisers. In its most recent iteration, `LaMa` allows for
automatic differentiation with the `RTMB` package which drastically
increases speed and accuracy even more.

The most important families of functions are

- the `forward` family that calculates the log-likelihood for various
  different models,

- the `tpm` family for calculating transition probability matrices,
  <!-- + `tpm()` for calculating a homogeneous transition probability matrix via the multinomial logistic link,  -->
  <!-- + `tpm_g()` for calculating general inhomogeneous transition probabilty matrices,  -->
  <!-- + `tpm_p()` for calculating transition matrices of periodically inhomogeneous HMMs, -->
  <!-- + `tpm_cont()` for calculating the transition probabilites of a continuous-time Markov chain, -->
  <!-- + `tpm_emb()` for calculating the embedded transition matrix of an HSMM, -->

- the `stationary` family to compute stationary and periodically
  stationary distributions

- as well as the `stateprobs` and `viterbi` families for local and
  global decoding.

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

- [Introduction to
  LaMa](https://janolefi.github.io/LaMa/articles/Intro_to_LaMa.html)
- [Inhomogeneous HMMs with covariate
  effects](https://janolefi.github.io/LaMa/articles/Inhomogeneous_HMMs.html)
- [Longitudinal
  data](https://janolefi.github.io/LaMa/articles/Longitudinal_data.html)
- [Periodic
  HMMs](https://janolefi.github.io/LaMa/articles/Periodic_HMM.html)
- [LaMa and
  RTMB](https://janolefi.github.io/LaMa/articles/LaMa_and_RTMB.html)
- [Penalised
  splines](https://janolefi.github.io/LaMa/articles/Penalised_splines.html)

Other latent Markov model classes:

- [State-space
  models](https://janolefi.github.io/LaMa/articles/State_space_models.html)
- [Continuous-time
  HMMs](https://janolefi.github.io/LaMa/articles/Continuous_time_HMMs.html)
- [Hidden semi-Markov
  models](https://janolefi.github.io/LaMa/articles/HSMMs.html)
- [Markov-modulated (marked) Poisson
  processes](https://janolefi.github.io/LaMa/articles/MMMPPs.html)

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
#>   0.415   0.019   0.471
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

hist(trex$step, prob = TRUE, bor = "white", breaks = 40, main = "", xlab = "step length")
curve(delta[1] * dgamma2(x, mu[1], sigma[1]), add = TRUE, lwd = 2, col = "orange", n=500)
curve(delta[2] * dgamma2(x, mu[2], sigma[2]), add = TRUE, lwd = 2, col = "deepskyblue", n=500)
legend("topright", col = c("orange", "deepskyblue"), lwd = 2, bty = "n", legend = c("state 1", "state 2"))
```

<img src="man/figures/README-visualization-1.png" width="75%" style="display: block; margin: auto;" />

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-mews2024build" class="csl-entry">

Mews, Sina, Jan-Ole Koslik, and Roland Langrock. 2025. “How to Build
Your Latent Markov Model: The Role of Time and Space.” *Statistical
Modelling* 0 (0). <https://doi.org/10.1177/1471082X251355681>.

</div>

<div id="ref-zucchini" class="csl-entry">

Zucchini, Walter, Iain L. MacDonald, and Roland Langrock. 2016. *Hidden
Markov Models for Time Series: An Introduction Using R*. Boca Raton:
Chapman & Hall/CRC.

</div>

</div>
