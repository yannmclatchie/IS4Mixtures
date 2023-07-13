library(ggplot2)
library(tidyverse)
library(posterior)
library(bayesflow)

# load utility functions
source("R/nnn/config.R")
source("R/nnn/utils.R")
source("R/nnn/proposals.R")
source("R/nnn/targets.R")
source("R/nnn/ggtheme.R")

# compute the normalising constant of the target distribution
Z <- integrate(linear_kernel, -Inf, Inf, weights = weights)
lZ <- log(Z$value)

# define some functional to estimate
functional <- function(x) sin(x) * cos(x)
functional <- Vectorize(functional)
I.quad <- integrate(function(x)
  functional(x) * target_dens(x, weights = weights, lZ = lZ),
  -Inf, Inf) 
I.quad$value # should be zero since functional is odd

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

# evaluate all proposals on linear target
N1_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = linear_target, 
                                  proposal = N1_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
N2_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = linear_target, 
                                  proposal = N2_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
N3_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = linear_target, 
                                  proposal = N3_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
complete_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = linear_target, 
                                  proposal = complete_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
aware_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = linear_target, 
                                  proposal = aware_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
res_df <- bind_rows(N1_proposal_res,
                    N2_proposal_res,
                    N3_proposal_res,
                    complete_proposal_res,
                    aware_proposal_res)
res_df

# plot the proposal performances on linear target
p_proposal_metrics <- res_df |> dplyr::select(proposal, bias, var, ess) |>
  group_by(proposal) |>
  reshape2::melt(id.vars = "proposal", variable.name = "metric") |>
  group_by(proposal, metric) |>
  summarise(mean_metric = mean(value),
            sd_metic = sd(value),
            lower_metric = quantile(value, prob = 0.05),
            upper_metric = quantile(value, prob = 0.95)) |>
  ggplot(aes(y = mean_metric, x = proposal)) +
  geom_pointrange(aes(ymin = lower_metric, 
                      ymax = upper_metric),
                  size = 0.4,
                  shape = 1) +
  facet_wrap(~metric, scales = "free") +
  ylab(NULL) +
  xlab(NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust=1))
p_proposal_metrics
#save_tikz_plot(p_proposal_metrics, width = 4, 
#               filename = "./tex/linear-nnn-proposals.tex")
