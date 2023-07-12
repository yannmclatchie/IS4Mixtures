library(ggplot2)
library(tidyverse)

# load utility functions
source("R/nnn/config.R")
source("R/nnn/utils.R")
source("R/nnn/proposals.R")
source("R/nnn/targets.R")
source("R/nnn/ggtheme.R")

# produce samples from the components
N1_samples <- rnorm(num_draws, mu_1, sigma_1)
N2_samples <- rnorm(num_draws, mu_2, sigma_2)
N3_samples <- rnorm(num_draws, mu_3, sigma_3)

# compute the normalising constant of the target distribution
Z <- integrate(linear_kernel, -Inf, Inf, weights = weights)
lZ <- log(Z$value)

# plot the (known) target and the component models
p <- ggplot() +
  geom_function(fun = linear_dens,
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
p

# define some functional to estimate
functional <- function(x) sin(x) * cos(x)
functional <- Vectorize(functional)
I.quad <- integrate(function(x)
  functional(x) * target_dens(x, weights = weights, lZ = lZ),
  -Inf, Inf) 
I.quad # should be zero since functional is odd

# investigate single-model proposals
proposal_name <- "N1"
draws <- N1_samples
lws_u <- (linear_log_kernel(draws, weights = weights) 
          - dnorm(draws, mu_1, sigma_1, log = TRUE))
ws_u <- exp(lws_u)
proposal_metrics(proposal_name, functional, draws, ws_u)

# investigate equally-weighted linear mixture proposal
# coming soon ...
