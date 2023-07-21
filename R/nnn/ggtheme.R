# define a ggplot theme
my_theme <- (theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.background = element_blank()))

# set ggplot theme
theme_set(my_theme)

# custom labeller
metric_names <- c(
  `var` = "Variance",
  `ess` = "ESS",
  `squared_error` = "Squared error",
  `est` = "Estimate"
)