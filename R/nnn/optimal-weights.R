# load utility functions
source("R/nnn/config.R")

optimal_weights <- function(m_k, v_k, alpha_k) {
  # compute some statistics
  l_k <- alpha_k / v_k
  m_ast <- sum(l_k * m_k) / sum(alpha_k)
  v_ast <- 1 / sum(l_k)
  
  # reparameterisations
  a_k <- (2 * v_k - v_ast) / (2 * v_ast * v_k)
  b_k <- (2 * v_k * m_ast - v_ast * m_k) / (2 * v_ast * v_k)
  c_k <- (2 * v_k * m_ast^2 - v_ast * m_k^2) / (2 * v_ast * v_k)

  # compute the estimator variances
  sigma2_k <- (exp(-c_k + b_k^2 / a_k) 
               * sqrt(2 * a_k * v_k) / v_ast 
               * (1 / (2 * a_k) + b_k^2 / a_k^2) - m_ast^2)
  
  # return the optimal weights
  W <- sum(1 / sigma2_k)
  return(1 / (W * sigma2_k))
}

## Scenario 1: all components are useful

# define the experiment 
m_k <- c(0, 1, 2)
v_k <- c(1, 1, 1)
alpha_k <- c(w1, w2, w3)

# compute the optimal weights
optimal_weights(m_k, v_k, alpha_k)

## Scenario 2: some of the components produce non-finite integrals

# define the experiment 
m_k <- c(mu_1, mu_2, mu_3)
v_k <- c(sigma_1^2, sigma_2^2, sigma_3^2)
alpha_k <- c(w1, w2, w3)

# compute the optimal weights
optimal_weights(m_k, v_k, alpha_k)
