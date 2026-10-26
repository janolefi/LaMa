## prototype single kernel: Gamma is N x N x nSlices, gidx[t] = slice for transition into t
proto_kernel <- function(delta, Gamma, allprobs, gidx, starts, ends) {
  "[<-" <- ADoverload("[<-")
  l <- 0
  for (i in seq_along(starts)) {
    s <- starts[i]; e <- ends[i]
    foo <- delta[i, , drop = FALSE] * allprobs[s, , drop = FALSE]
    sumfoo <- sum(foo); phi <- foo / sumfoo; l <- l + log(sumfoo)
    for (t in seq_len(e - s) + s) {
      foo <- (phi %*% Gamma[, , gidx[t]]) * allprobs[t, , drop = FALSE]
      sumfoo <- sum(foo); phi <- foo / sumfoo; l <- l + log(sumfoo)
    }
  }
  l
}
## row-max shift wrapper (logspace input)
proto_shift <- function(delta, Gamma, logP, gidx, starts, ends) {
  m <- apply(logP, 1, max)
  proto_kernel(delta, Gamma, exp(logP - m), gidx, starts, ends) + sum(m)
}
