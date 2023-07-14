# single-model proposals
N1_proposal <- list(name = "N1",
                    sampling_weights = c(rep(1, length(N1_samples)),
                                         rep(0, length(N2_samples)),
                                         rep(0, length(N3_samples))),
                    log_kernel = function(x, mu_1, sigma_1, ...){
                      -1 / (2 * sigma_1^2) * (x - mu_1)^2},
                    log_dens = function(x, mu_1, sigma_1, ...) {
                      dnorm(x, mu_1, sigma_1, log = TRUE)})

N2_proposal <- list(name = "N2",
                    sampling_weights = c(rep(0, length(N1_samples)),
                                         rep(1, length(N2_samples)),
                                         rep(0, length(N3_samples))),
                    log_kernel = function(x, mu_2, sigma_2, ...){
                      -1 / (2 * sigma_2^2) * (x - mu_2)^2},
                    log_dens = function(x, mu_2, sigma_2, ...) {
                      dnorm(x, mu_2, sigma_2, log = TRUE)})

N3_proposal <- list(name = "N3",
                    sampling_weights = c(rep(0, length(N1_samples)),
                                         rep(0, length(N2_samples)),
                                         rep(1, length(N3_samples))),
                    log_kernel = function(x, mu_3, sigma_3, ...){
                      -1 / (2 * sigma_3^2) * (x - mu_3)^2},
                    log_dens = function(x, mu_3, sigma_3, ...) {
                      dnorm(x, mu_3, sigma_3, log = TRUE)})

# equally-weighted linear combination
complete_log_kernel <- function(x, ...) {
  log(1 / 3 * dnorm(x, mu_1, sigma_1) 
      + 1 / 3 * dnorm(x, mu_2, sigma_2) 
      + 1 / 3 * dnorm(x, mu_3, sigma_3))
}
complete_kernel <- function(x, ...) exp(complete_log_kernel(x, ...))
complete_Z <- integrate(complete_kernel, -Inf, Inf)
complete_lZ <- log(complete_Z$value)
complete_log_dens <- function(x, ...) complete_log_kernel(x, ...) - complete_lZ
complete_proposal <- list(name = "complete",
                          sampling_weights = c(rep(1, length(N1_samples)),
                                               rep(1, length(N2_samples)),
                                               rep(1, length(N3_samples))),
                          log_kernel = complete_log_kernel,
                          log_dens = complete_log_dens)

# structure-aware combination
aware_log_kernel <- function(x, ...) {
  log(weights$w1 * dnorm(x, mu_1, sigma_1) 
      + weights$w2 * dnorm(x, mu_2, sigma_2) 
      + weights$w3 * dnorm(x, mu_3, sigma_3))
}
aware_kernel <- function(x, ...) exp(aware_log_kernel(x, ...))
aware_Z <- integrate(aware_kernel, -Inf, Inf, weights = weights)
aware_lZ <- log(aware_Z$value)
aware_log_dens <- function(x, ...) aware_log_kernel(x, ...) - aware_lZ
aware_proposal <- list(name = "aware",
                        sampling_weights = c(rep(weights$w1, length(N1_samples)),
                                             rep(weights$w2, length(N2_samples)),
                                             rep(weights$w3, length(N3_samples))),
                        log_kernel = aware_log_kernel,
                        log_dens = aware_log_dens)
