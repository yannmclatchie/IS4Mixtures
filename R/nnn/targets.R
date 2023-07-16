# linear stacking
linear_kernel <- function(x, weights, ...) {
  (weights$w1 * dnorm(x, mu_1, sigma_1) 
   + weights$w2 * dnorm(x, mu_2, sigma_2)
   + weights$w3 * dnorm(x, mu_3, sigma_3))
}
linear_log_kernel <- function(x, weights, ...) {
  log(linear_kernel(x, weights))
}
linear_Z <- integrate(linear_kernel, -Inf, Inf, weights = weights)
linear_lZ <- log(linear_Z$value)
linear_log_dens <- function(x, weights, ...) {
  linear_log_kernel(x = x, weights = weights) - linear_lZ
}
linear_dens <- function(x, weights, ...) exp(linear_log_dens(x, weights))
linear_target <- list(kernel = linear_kernel,
                      log_kernel = linear_log_kernel,
                      log_dens = linear_log_dens,
                      dens = linear_dens)

# log-linear stacking
locking_log_kernel <- function(x, weights, ...) {
  (weights$w1 * dnorm(x, mu_1, sigma_1, log = TRUE) 
   + weights$w2 * dnorm(x, mu_2, sigma_2, log = TRUE)
   + weights$w3 * dnorm(x, mu_3, sigma_3, log = TRUE))
}
locking_kernel <- function(x, weights, ...) {
  exp(locking_log_kernel(x, weights))
}
locking_Z <- integrate(locking_kernel, -Inf, Inf, weights = weights)
locking_lZ <- log(locking_Z$value)
locking_log_dens <- function(x, weights, ...) {
  locking_log_kernel(x = x, weights = weights) - locking_lZ
}
locking_dens <- function(x, weights, ...) exp(locking_log_dens(x, weights))
locking_target <- list(kernel = locking_kernel,
                       log_kernel = locking_log_kernel,
                       log_dens = locking_log_dens,
                       dens = locking_dens)

# quantum super-position
quacking_log_kernel <- function(x, weights, w0, betas) {
  # Equation 8 of Yao et al. (2022)
  w0 * log(betas$b1 * dnorm(x, mu_1, sigma_1) 
           + betas$b2 * dnorm(x, mu_2, sigma_2)
           + betas$b3 * dnorm(x, mu_3, sigma_3)) + locking_log_kernel(x, weights)
}
quacking_kernel <- function(x, weights, w0, betas) {
  exp(quacking_log_kernel(x, weights, w0, betas))
}
quacking_Z <- integrate(quacking_kernel, -Inf, Inf, weights = weights,
                        w0 = w0, betas = betas)
quacking_lZ <- log(quacking_Z$value)
quacking_log_dens <- function(x, weights, w0, betas) {
  quacking_log_kernel(x = x, weights = weights, 
                      w0 = w0, betas = betas) - quacking_lZ
}
quacking_dens <- function(x, weights, w0, betas) {
  exp(quacking_log_dens(x, weights, w0, betas))
}
quacking_target <- list(kernel = quacking_kernel,
                        log_kernel = quacking_log_kernel,
                        log_dens = quacking_log_dens,
                        dens = quacking_dens)
