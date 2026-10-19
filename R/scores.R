# Number of observations in all forward() calls of nll (counting mode),
# found by one evaluation of nll at par without AD
count_obs <- function(nll, par){
  .lama_env$count <- 0
  on.exit(.lama_env$count <- NULL)
  nll(par)
  .lama_env$count
}
