# indicator function to test coverage
is_in <- function(x, l, u) {
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  return(result)
}

# compute the ESS of a set of importance weights
ess_IS <- function(ws) {
  S <- sum(ws)
  return(S^2 / sum(ws^2))
}

# compute the estimate and variance of a function by IS
proposal_metrics <- function(proposal_name, functional, draws, 
                             ws_u, level = .95) {
  # normalise the (un-normalised) weights
  ws <- ws_u / sum(ws_u)
  # compute the functional at the sampled draws
  fs <- functional(draws)
  # compute the ESS, mean, and variance of the IS estimate
  ess <- ess_IS(ws)
  mean <-  sum(fs * ws)
  var <- sum(ws^2 * (fs - mean)^2)
  # compute the quartiles corresponding to the level
  D <- qnorm(p = (1 + level) / 2)
  d <- D * sqrt(var)
  # return the metrics
  list(proposal = proposal_name,
       est = mean,
       var = var,
       lwr = mean - d,
       upr = mean + d,
       ess = ess)
}
