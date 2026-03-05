# create environment to store information of distributions and parameters used in allprobs calculations
# forward() and forward_g() then access this environment and REPORT() the contents 
# this way we can capture all the parameters used in the observation distributions
# .obs_info <- new.env(parent = emptyenv())

# dist_reporter <- function(fun, family = NULL) {
#   if(is.null(family)) {
#     family <- substring(deparse(substitute(fun)), 2)
#   } 
#   
#   function(...) {
#     mc <- match.call(fun, expand.dots = FALSE)
#     stream <- deparse(mc[[2]])
#     
#     # match call to formal arguments
#     matched_call <- match.call(fun, call = mc)
#     matched_list <- as.list(matched_call)[-1]           # drop function name
#     matched_list <- matched_list[names(matched_list) != "x"]  # drop data argument
#     
#     # evaluate data and parameters
#     eval_args <- lapply(matched_list, eval, parent.frame())
#     obs_data <- eval(mc[[2]], parent.frame())  # get the actual observed data
#     
#     # initialize storage if first call
#     if (is.null(.obs_info$calls[[stream]])) {
#       .obs_info$calls[[stream]] <- list(
#         family = family,
#         obs    = obs_data,                   # store data once
#         pars   = lapply(eval_args, function(x) list())  # list per parameter
#       )
#     }
#     
#     # append each parameter value as a new element in the list
#     for (nm in names(eval_args)) {
#       .obs_info$calls[[stream]]$pars[[nm]][[length(.obs_info$calls[[stream]]$pars[[nm]]) + 1]] <- eval_args[[nm]]
#     }
#     
#     # call the original function
#     do.call(fun, list(...))
#   }
# }