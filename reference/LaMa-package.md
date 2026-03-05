# LaMa: Fast Numerical Maximum Likelihood Estimation for Latent Markov Models

A variety of latent Markov models, including hidden Markov models,
hidden semi-Markov models, state-space models and continuous-time
variants can be formulated and estimated within the same framework via
directly maximising the likelihood function using the so-called forward
algorithm. Applied researchers often need custom models that standard
software does not easily support. Writing tailored 'R' code offers
flexibility but suffers from slow estimation speeds. We address these
issues by providing easy-to-use functions (written in 'C++' for speed)
for common tasks like the forward algorithm. These functions can be
combined into custom models in a Lego-type approach, offering up to
10-20 times faster estimation via standard numerical optimisers. To aid
in building fully custom likelihood functions, several vignettes are
included that show how to simulate data from and estimate all the above
model classes.

## See also

Useful links:

- <https://janolefi.github.io/LaMa/>

## Author

**Maintainer**: Jan-Ole Fischer <jan-ole.fischer@mailbox.org>
([ORCID](https://orcid.org/0009-0004-1556-9053))
