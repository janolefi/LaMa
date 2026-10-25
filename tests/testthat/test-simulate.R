test_that("sim_states(): initial distribution, transitions and tracks", {
  set.seed(1)
  Gamma <- tpm(c(-2, -1.5)); delta <- c(0.9, 0.1)
  inp <- prepare_decoding(delta, Gamma, matrix(1, 20000, 2), rep(1:200, each = 100))
  Z <- sim_states(inp)
  expect_equal(mean(Z[inp$start] == 1), 0.9, tolerance = 0.05)
  inner <- setdiff(2:20000, inp$start)
  tr <- table(Z[inner - 1], Z[inner]) / as.vector(table(Z[inner - 1]))
  expect_equal(unname(tr), unname(Gamma), tolerance = 0.02)
})

test_that("sim_states(): time-varying Gamma uses the slice into t", {
  set.seed(2)
  n <- 40000
  G <- array(0, c(2, 2, n))
  for(t in 1:n) G[, , t] <- if(t %% 2) tpm(c(-4, -4)) else tpm(c(4, 4)) # odd: stay, even: switch
  Z <- sim_states(prepare_decoding(c(0.5, 0.5), G, matrix(1, n, 2)))
  expect_gt(mean(Z[seq(3, n, 2)] == Z[seq(2, n - 1, 2)]), 0.95)
  expect_gt(mean(Z[seq(2, n, 2)] != Z[seq(1, n - 1, 2)]), 0.95)
})

test_that("sim_stream(): per state, per observation, missing values, r-function lookup", {
  set.seed(3)
  n <- 20000
  Z <- rep(1:2, length.out = n)
  x <- rep(1, n); x[7] <- NA
  st <- list(x = x, dens = dnorm, name = "dnorm", pars = list(c(0, 10), 1), nStates = 2)
  y <- sim_stream(st, Z)
  expect_true(is.na(y[7]) && sum(is.na(y)) == 1)
  expect_equal(c(mean(y[Z == 1], na.rm = TRUE), mean(y[Z == 2], na.rm = TRUE)), c(0, 10), tolerance = 0.05)
  # parameters per observation and state, and a different namespace
  st$pars[[1]] <- cbind(seq(0, 1, length = n), seq(10, 11, length = n))
  y <- sim_stream(st, Z)
  expect_equal(mean(y[Z == 2 & seq_len(n) > n / 2], na.rm = TRUE), 10.75, tolerance = 0.05)
  st$dens <- dgamma2; st$name <- "LaMa::dgamma2" # from another namespace, name with prefix
  st$pars <- list(c(1, 5), c(0.5, 0.5))
  expect_equal(mean(sim_stream(st, Z)[Z == 2], na.rm = TRUE), 5, tolerance = 0.05)
  # no counterpart: error, unless given
  zig <- zero_inflate(dnorm)
  st$dens <- zig; st$name <- "zig"; st$pars <- list(0, c(0, 10), 1)
  expect_error(sim_stream(st, Z), "rdens")
  rzig <- function(n, zeroprob, mean, sd) rnorm(n, mean, sd)
  expect_equal(mean(sim_stream(st, Z, rzig)[Z == 2], na.rm = TRUE), 10, tolerance = 0.05)
})
