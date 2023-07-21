# set the seed for reproducibility
SEED <- 1234
set.seed(SEED)

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
w1 <- 0.85
w2 <- 0.05
w3 <- 0.1
w1 + w2 + w3

# build weights list
weights <- list(w1 = w1, w2 = w2, w3 = w3)

# define quacking-specific parameters
w0 <- 0.5
b1 <- 0.6
b2 <- 0.0
b3 <- 0.4
betas <- list(b1 = b1, b2 = b2, b3 = b3)

# define a function to estimate
functional <- function(x) x
functional <- Vectorize(functional)
