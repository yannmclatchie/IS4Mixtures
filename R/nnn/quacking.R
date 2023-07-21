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

# estimate the functional under the quacked kernel
I.quad <- integrate(function(x){
  functional(x) * quacking_dens(x, weights = weights,
                                w0 = w0, betas = betas)},
  -Inf, Inf) 
I.quad$value

# plot the (known) target and the component models
p_quacking_target <- ggplot() +
  geom_function(fun = quacking_dens,
                args = list(weights = weights,
                            w0 = w0,
                            betas = betas)) + 
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
  ) +
  xlab(NULL) + 
  ylab(NULL)
p_quacking_target
#save_tikz_plot(p_quacking_target, width = 3, 
#               filename = "./tex/quacking-nnn-target.tex")

# evaluate all proposals on locking target
N1_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = quacking_target, 
                                  proposal = N1_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
N2_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = quacking_target, 
                                  proposal = N2_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
N3_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = quacking_target, 
                                  proposal = N3_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
complete_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = quacking_target, 
                                  proposal = complete_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
aware_proposal_res <- 1:num_iters |>
  map(\(rep_id) rep_eval_proposal(rep_id, target = quacking_target, 
                                  proposal = aware_proposal, 
                                  I.quad = I.quad$value, 
                                  functional = functional)) |> 
  bind_rows()
quacking_res_df <- bind_rows(N1_proposal_res,
                            N2_proposal_res,
                            N3_proposal_res,
                            complete_proposal_res,
                            aware_proposal_res)
quacking_res_df

# plot the proposal performances on linear target
p_quacking_metrics <- quacking_res_df |> dplyr::select(proposal, snis, error, var, ess) |>
  mutate(squared_error = error^2) |>
  dplyr::select(-error) |>
  #dplyr::select(-squared_error) |> # uncomment this line to remove squared error
  reshape2::melt(id.vars = c("proposal", "snis"), variable.name = "metric") |>
  group_by(proposal, snis, metric) |>
  summarise(mean_metric = mean(value),
            sd_metic = sd(value),
            lower_metric = quantile(value, prob = 0.05),
            upper_metric = quantile(value, prob = 0.95)) |>
  ggplot(aes(y = mean_metric, x = proposal, colour = snis)) +
  geom_pointrange(aes(ymin = lower_metric, 
                      ymax = upper_metric),
                  size = 0.4,
                  shape = 1,
                  position = position_dodge(width = 0.4)) +
  facet_wrap(~metric, scales = "free", labeller = as_labeller(metric_names)) +
  ylab(NULL) +
  xlab(NULL) +
  scale_y_log10() +
  scale_colour_manual(values = c("black", "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        legend.position = "none")
p_quacking_metrics
#save_tikz_plot(p_quacking_metrics, width = 5, 
#               filename = "./tex/quacking-nnn-proposals.tex")

