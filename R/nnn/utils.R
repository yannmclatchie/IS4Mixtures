# evaluate a given proposal for a target
eval_proposal <- function(target, proposal, all_draws, I.quad, functional) {
  # unpack inputs
  target_kernel <- target$kernel
  target_log_kernel <- target$log_kernel
  target_dens <- target$dens
  proposal_name <- proposal$name
  proposal_sampling_weights <- proposal$sampling_weights
  proposal_log_kernel <- proposal$log_kernel
  
  # resample from the collection of all draws
  draws_rs <- resample_draws(all_draws, 
                             weights = proposal_sampling_weights,
                             method = "simple_no_replace",
                             ndraws = num_draws)
  
  # compute the importance ratios from the resampled draws
  lws_u <- (target_log_kernel(draws_rs$x, weights = weights) 
            - proposal_log_kernel(draws_rs$x, 
                                  mu_1 = mu_1, sigma_1 = sigma_1, 
                                  mu_2 = mu_2, sigma_2 = sigma_2, 
                                  mu_3 = mu_3, sigma_3 = sigma_3))
  ws_u <- exp(lws_u)
  
  # output the proposal metrics
  proposal_metrics(proposal_name, functional, I.quad, draws, ws_u)
}

# simulate the evalution for one iteration 
rep_eval_proposal <- function(rep_id, target, proposal, 
                              all_draws, I.quad, functional) {
  # produce samples from the components
  N1_samples <- rnorm(num_draws, mu_1, sigma_1)
  N2_samples <- rnorm(num_draws, mu_2, sigma_2)
  N3_samples <- rnorm(num_draws, mu_3, sigma_3)
  
  # concatenate draws
  all_draws <- as_draws_df(list(x = c(N1_samples, N2_samples, N3_samples)))
  
  # evaluate 
  out <- eval_proposal(target, proposal, all_draws, I.quad, functional)
  out$rep = rep_id
  return(out)
}

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
proposal_metrics <- function(proposal_name, functional, I.quad,
                             draws, ws_u, level = .95) {
  # normalise the (un-normalised) weights
  ws <- ws_u / sum(ws_u)
  # compute the functional at the sampled draws
  fs <- functional(draws)
  # compute the ESS, mean, and variance of the IS estimate
  ess <- ess_IS(ws)
  mean <-  sum(fs * ws)
  var <- sum(ws^2 * (fs - mean)^2)
  error <- I.quad - mean
  # compute the quartiles corresponding to the level
  D <- qnorm(p = (1 + level) / 2)
  d <- D * sqrt(var)
  # return the metrics
  list(proposal = proposal_name,
       est = mean,
       var = var,
       error = error,
       lwr = mean - d,
       upr = mean + d,
       ess = ess)
}

