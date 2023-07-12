# set the seed
SEED <- 1234

# define the number of draws (use default Stan)
num_draws <- 4e3
num_iters <- 500

# define the components
mu_1 <- 20
mu_2 <- 2
mu_3 <- 0
sigma_1 <- 5
sigma_2 <- 3
sigma_3 <- 1

# define the combination weights
w1 <- 0.85#0.65
w2 <- 0.05
w3 <- 0.1#0.3
w1 + w2 + w3

# build weights list
weights <- list(w1 = w1, w2 = w2, w3 = w3)
