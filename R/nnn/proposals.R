# single-model proposals
N1_proposal <- list(name = "N1",
                    sampling_weights = c(rep(1, length(N1_samples)),
                                         rep(0, length(N2_samples)),
                                         rep(0, length(N3_samples))),
                    log_kernel = function(x, mu_1, sigma_1, ...) {
                    dnorm(x, mu_1, sigma_1, log = TRUE)})

N2_proposal <- list(name = "N2",
                    sampling_weights = c(rep(0, length(N1_samples)),
                                         rep(1, length(N2_samples)),
                                         rep(0, length(N3_samples))),
                    log_kernel = function(x, mu_2, sigma_2, ...) {
                      dnorm(x, mu_2, sigma_2, log = TRUE)})

N3_proposal <- list(name = "N3",
                    sampling_weights = c(rep(0, length(N1_samples)),
                                         rep(0, length(N2_samples)),
                                         rep(1, length(N3_samples))),
                    log_kernel = function(x, mu_3, sigma_3, ...) {
                      dnorm(x, mu_3, sigma_3, log = TRUE)})

# equally-weighted linear combination
complete_log_kernel <- function(x, ...) {
  log(dnorm(x, mu_1, sigma_1) 
      + dnorm(x, mu_2, sigma_2) 
      + dnorm(x, mu_3, sigma_3))
}
complete_proposal <- list(name = "complete",
                          sampling_weights = c(rep(1, length(N1_samples)),
                                               rep(1, length(N2_samples)),
                                               rep(1, length(N3_samples))),
                          log_kernel = complete_log_kernel)

# structure-aware combination
aware_log_kernel <- function(x, ...) {
  log(weights$w1 * dnorm(x, mu_1, sigma_1) 
      + weights$w2 * dnorm(x, mu_2, sigma_2) 
      + weights$w3 * dnorm(x, mu_3, sigma_3))
}
aware_proposal <- list(name = "aware",
                        sampling_weights = c(rep(weights$w1, length(N1_samples)),
                                             rep(weights$w2, length(N2_samples)),
                                             rep(weights$w3, length(N3_samples))),
                        log_kernel = aware_log_kernel)