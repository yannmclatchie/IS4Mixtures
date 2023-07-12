# equally-weighted linear combination
equal_linear_dens <- function(x) {
  dnorm(x, mu_1, sigma_1) + dnorm(x, mu_2, sigma_2) + dnorm(x, mu_3, sigma_3)
}