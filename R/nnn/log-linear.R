library(ggplot2)
library(tidyverse)

source("R/nnn/config.R")
source("R/nnn/proposals.R")
source("R/nnn/targets.R")
source("R/nnn/ggtheme.R")

# produce samples from the components
N1_samples <- rnorm(num_draws, mu_1, sigma_1)
N2_samples <- rnorm(num_draws, mu_2, sigma_2)
N3_samples <- rnorm(num_draws, mu_3, sigma_3)

# plot the (known) target and the component models
Z <- integrate(locking_kernel,-Inf, Inf, weights = weights)
lZ <- log(Z$value)
ggplot() +
  geom_function(fun = locking_dens,
                args = list(weights = weights, lZ = lZ)) + 
  geom_density(
    data = data.frame(x = N1_samples), aes(x),
    fill = "grey20", colour = NA, alpha = 0.1
  ) +
  geom_density(
    data = data.frame(x = N2_samples), aes(x),
    fill = "grey20", colour = NA, alpha = 0.1
  ) +
  geom_density(
    data = data.frame(x = N3_samples), aes(x),
    fill = "grey20", colour = NA, alpha = 0.1
  )

