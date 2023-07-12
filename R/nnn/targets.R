# linear stacking
linear_kernel <- function(x, weights) {
  (weights$w1 * dnorm(x, mu_1, sigma_1) 
   + weights$w2 * dnorm(x, mu_2, sigma_2)
   + weights$w3 * dnorm(x, mu_3, sigma_3))
}
linear_log_kernel <- function(x, weights) {
  log(linear_kernel(x, weights))
}
linear_dens <- function(x, weights, lZ) {
  ldens <- linear_log_kernel(x = x, weights = weights) - lZ
  return(exp(ldens))
}

# log-linear stacking
locking_log_kernel <- function(x, weights) {
  (weights$w1 * dnorm(x, mu_1, sigma_1, log = TRUE) 
   + weights$w2 * dnorm(x, mu_2, sigma_2, log = TRUE)
   + weights$w3 * dnorm(x, mu_3, sigma_3, log = TRUE))
}
locking_kernel <- function(x, weights) {
  exp(locking_log_kernel(x, weights))
}
locking_dens <- function(x, weights, lZ) {
  ldens <- locking_log_kernel(x = x, weights = weights) - lZ
  return(exp(ldens))
}

# quantum super-position
# (coming soon ... maybe)

