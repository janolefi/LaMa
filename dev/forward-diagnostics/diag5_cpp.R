source(file.path(Sys.getenv("SP"), "helpers.R"))
Rcpp::sourceCpp(file.path(Sys.getenv("SP"), "proto.cpp"))
G1 <- matrix(c(.9, .1, .2, .8), 2, byrow = TRUE); d0 <- c(2/3, 1/3)
n <- 1e6; lp <- sim_logP(n, N = 3); P <- exp(lp)
G3 <- tpm(rep(-2, 6)); D <- matrix(rep(1/3, 3), 1)
Gt <- array(G3, c(3, 3, n)); Gp <- array(G3, c(3, 3, 24)); tod <- rep(1:24, length.out = n)
tm <- function(e) unname(system.time(for (i in 1:5) force(e))[3] / 5)
cat("n = 1e6, N = 3, per call seconds\n")
cat("homog   old:", tm(LaMa:::forward_cpp_h(P, D[1, ], G3)), " new:", tm(forward_cpp_proto(P, D, array(G3, c(3,3,1)), rep(0L, n), 0L)), "\n")
cat("inhom   old:", tm(LaMa:::forward_cpp_g(P, D[1, ], Gt[, , -1])), " new:", tm(forward_cpp_proto(P, D, Gt, 0:(n - 1), 0L)), "\n")
cat("periodic old (expand Gamma[,,tod] + g):", tm(LaMa:::forward_cpp_g(P, D[1, ], Gp[, , tod][, , -1])),
    " new:", tm(forward_cpp_proto(P, D, Gp, tod - 1L, 0L)), "\n")
cat("full forward() R wrapper, homog:", tm(forward(D[1, ], G3, P, ad = FALSE, report = FALSE)), "\n")
cat("full forward_p() R wrapper:", tm(forward_p(D[1, ], Gp, P, tod, ad = FALSE, report = FALSE)), "\n")
ID <- rep(1:100, each = n / 100); st <- which(!duplicated(ID)) - 1L
cat("tracks  old:", tm(LaMa:::forward_cpp_h_tracks(P, D[rep(1, 100), ], array(G3, c(3, 3, 100)), calc_trackInd(ID))),
    " new:", tm(forward_cpp_proto(P, D[rep(1, 100), ], array(G3, c(3, 3, 1)), rep(0L, n), st)), "\n")
cat("row max via apply:", tm(apply(lp, 1, max)), "  via pmax loop:", tm({m <- lp[, 1]; for (j in 2:3) m <- pmax(m, lp[, j]); m}), "\n")
cat("equal:", all.equal(LaMa:::forward_cpp_g(P, D[1, ], Gp[, , tod][, , -1]), forward_cpp_proto(P, D, Gp, tod - 1L, 0L)), "\n")
