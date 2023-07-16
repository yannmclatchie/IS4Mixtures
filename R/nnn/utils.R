# evaluate a given proposal for a target
eval_proposal <- function(rep_id, target, proposal, all_draws, 
                          I.quad, functional, self_norm) {
  # unpack inputs
  target_kernel <- target$kernel
  target_log_kernel <- target$log_kernel
  target_log_dens <- target$log_dens
  target_dens <- target$dens
  proposal_name <- proposal$name
  proposal_sampling_weights <- proposal$sampling_weights
  proposal_log_kernel <- proposal$log_kernel
  proposal_log_dens <- proposal$log_dens
  
  # resample from the collection of all draws
  draws_rs <- resample_draws(all_draws, 
                             weights = proposal_sampling_weights,
                             method = "simple_no_replace",
                             ndraws = num_draws)
  
  # compute the normalising constant of the target distribution
  Z <- integrate(target_kernel, -Inf, Inf, weights = weights, w0 = w0, betas = betas)
  lZ <- log(Z$value)
  
  # compute the importance ratios from the resampled draws
  lws_u <- (target_log_dens(draws_rs$x, weights = weights, 
                            w0 = w0, betas = betas) 
            - proposal_log_dens(draws_rs$x, 
                                mu_1 = mu_1, sigma_1 = sigma_1, 
                                mu_2 = mu_2, sigma_2 = sigma_2, 
                                mu_3 = mu_3, sigma_3 = sigma_3))
  ws_u <- exp(lws_u)
  
  # compute the proposal metrics
  res <- proposal_metrics(proposal_name, functional, I.quad, draws, ws_u)
  res$rep <- rep_id
  
  if (self_norm) {
    # and the SNIS importance ratios
    lws_sn_u <- (target_log_kernel(draws_rs$x, weights = weights,
                                   w0 = w0, betas = betas) 
                 - proposal_log_kernel(draws_rs$x, 
                                       mu_1 = mu_1, sigma_1 = sigma_1, 
                                       mu_2 = mu_2, sigma_2 = sigma_2, 
                                       mu_3 = mu_3, sigma_3 = sigma_3))
    ws_sn_u <- exp(lws_sn_u)
    
    # append SNIS metrics 
    snis_res <- proposal_metrics_sn(proposal_name, functional, I.quad, draws, ws_u)
    snis_res$rep <- rep_id
    
    res <- list(snis_res, res)
  }
  
  # return proposal metrics
  return(res)
}

# simulate the evalution for one iteration 
rep_eval_proposal <- function(rep_id, target, proposal, 
                              all_draws, I.quad, functional,
                              self_norm = TRUE) {
  # produce samples from the components
  # use the rep_id as the seed for reproducibility
  N1_samples <- withr::with_seed(rep_id, rnorm(num_draws, mu_1, sigma_1))
  N2_samples <- withr::with_seed(rep_id, rnorm(num_draws, mu_2, sigma_2))
  N3_samples <- withr::with_seed(rep_id, rnorm(num_draws, mu_3, sigma_3))
  
  # concatenate draws
  all_draws <- as_draws_df(list(x = c(N1_samples, N2_samples, N3_samples)))
  
  # evaluate 
  out <- eval_proposal(rep_id, target, proposal, all_draws, I.quad, 
                       functional, self_norm)
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

# compute the estimate and variance of a function by IS and SNIS
proposal_metrics_sn <- function(proposal_name, functional, I.quad,
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
       snis = "true",
       est = mean,
       var = var,
       error = error,
       lwr = mean - d,
       upr = mean + d,
       ess = ess)
}
proposal_metrics <- function(proposal_name, functional, I.quad,
                             draws, ws_u, level = .95) {
  # normalise the (un-normalised) weights
  ws <- ws_u / sum(ws_u)
  # compute the functional at the sampled draws
  fs <- functional(draws)
  # compute the ESS, mean, and variance of the IS estimate
  ess <- ess_IS(ws)
  mean <-  mean(fs * ws)
  var <- mean(ws^2 * (fs - mean)^2)
  error <- I.quad - mean
  # compute the quartiles corresponding to the level
  D <- qnorm(p = (1 + level) / 2)
  d <- D * sqrt(var)
  # return the metrics
  list(proposal = proposal_name,
       snis = "false",
       est = mean,
       var = var,
       error = error,
       lwr = mean - d,
       upr = mean + d,
       ess = ess)
}

