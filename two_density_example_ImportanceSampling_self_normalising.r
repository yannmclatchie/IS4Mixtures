f1_u <- function(x, log = FALSE) {
  ans <- -x ^ 2 / (2 * sigmasq)
  if (!log)
    ans <- exp(ans)
  return(ans)
}
f1_u <- Vectorize(f1_u)
rf1 <- function(n) {
  rnorm(n = n, sd = sqrt(sigmasq))
}
f2_u <- function(x, log = FALSE) {
  ans <- -(nu + 1) / 2 * log1p(x ^ 2 / nu)
  if (!log)
    ans <- exp(ans)
  return(ans)
}
f2_u <- Vectorize(f2_u)

rf2 <- function(n) {
  rt(n = n, df = nu)
}

pi_kernel <- function(x, weight, log = FALSE) {
  ans <- weight * f1_u(x, log = TRUE) +
    (1 - weight) * f2_u(x, log = TRUE)
  if (!log)
    ans <- exp(ans)
  return(ans)
}
pi_kernel <- Vectorize(pi_kernel)

pi <- function(x, weight, lZ) {
  ldens <- pi_kernel(x = x, weight = weight, log = TRUE) - lZ
  return(exp(ldens))
}
pi <- Vectorize(pi)

#TODO: implement log-weights and logSumExp-based estimators
w1_u <- function(x, weight) {
  logISw <- (1 - weight) * (f2_u(x, log = TRUE) -
                              f1_u(x, log = TRUE))
  return(exp(logISw))
}
w1_u <- Vectorize(w1_u)

w2_u <- function(x, weight) {
  logISw <- weight * (f1_u(x, log = TRUE) -
                        f2_u(x, log = TRUE))
  return(exp(logISw))
}
w2_u <- Vectorize(w2_u)

ess_IS <- function(ws) {
  S <- sum(ws)
  return(S ^ 2 / sum(ws ^ 2))
}

compute_both_estimators_SN <- function(M, weight, level = .95) {
  X.samp <- rf1(M)
  f1s <- functional(X.samp)
  W1s_u <- w1_u(X.samp, weight = weight)
  S1 <- sum(W1s_u)
  W1s <- W1s_u / S1
  ess.1 <- ess_IS(W1s)
  mu1 <-  sum(f1s * W1s)
  v1 <- sum(W1s ^ 2 * (f1s - mu1) ^ 2)
  #
  Y.samp <- rf2(M)
  f2s <- functional(Y.samp)
  W2s_u <- w2_u(Y.samp, weight = weight)
  S2 <- sum(W2s_u)
  W2s <- W2s_u / S2
  ess.2 <- ess_IS(W2s)
  mu2  <- sum(f2s * W2s)
  v2 <- sum(W2s ^ 2 * (f2s - mu2) ^ 2)
  #
  D <- qnorm(p = (1 + level) / 2)
  d1 <- D * sqrt(v1)
  d2 <- D * sqrt(v2)
  out1 <- tibble::tibble(
    est = mu1,
    var = v1,
    lwr = mu1 - d1,
    upr = mu1 + d1,
    ess = ess.1,
    proposal = "f[1]"
  )
  out2 <- tibble::tibble(
    est = mu2,
    var = v2,
    lwr = mu2 - d2,
    upr = mu2 + d2,
    ess = ess.2,
    proposal = "f[2]"
  )
  return(rbind(out1, out2))
}

is_in <- function(x, l, u) {
  below <- x >= l
  above <- x <= u
  result <- as.logical(below * above)
  return(result)
}

###########
sigmasq <- 1 ^ 2
nu <- 1
theta <- .01

Z <- integrate(pi_kernel,-Inf, Inf, weight = theta)
logZtheta <- log(Z$value)

minx <- -5
maxx <- 5

mY <-
  max(c(f1_u(0), f2_u(0), pi(0, weight = theta, lZ = logZtheta)))

curve(f2_u,
      minx,
      maxx,
      col = 2,
      lwd = 3,
      ylim = c(0, mY))
curve(f1_u, minx, maxx, lwd = 3, add = TRUE)
curve(
  pi(x, weight = theta, lZ = logZtheta),
  minx,
  maxx,
  col = "lightblue",
  lwd = 3,
  add = TRUE
)
legend(
  x = "topright",
  legend = c(expression(f[1]),
             expression(f[2]),
             expression(pi)),
  lwd = 2,
  col = c(1, 2, "lightblue"),
  bty = 'n'
)

curve(
  w1_u(x, weight = theta),
  minx,
  maxx,
  lwd = 3,
  log = "y",
  ylab = "Unnormalised importance weight"
)
curve(
  w2_u(x, weight = theta),
  minx,
  maxx,
  log = "y",
  col = 2,
  lty = 2,
  lwd = 3,
  add = TRUE
)
legend(
  x = "top",
  legend = c(expression(w[1]), expression(w[2])),
  lwd = 2,
  lty = c(1, 2),
  col = c(1, 2),
  bty = 'n'
)

###
## checking the algebra
functional <- function(x)
  x# sin(x)^2 + cos(x)
functional <- Vectorize(functional)

curve(
  functional,
  minx,
  maxx,
  lwd = 3,
  xlab = expression(X),
  ylab = expression(psi(x))
)

I.quad <- integrate(function(x)
  functional(x) * pi(x, weight = theta, lZ = logZtheta),-Inf,
  Inf) ## for sin(x) we know it should be zero (even function over the line)

I.quad

Nrep <- 500
Ndraws <- 400 # S, number of 'posterior' samples
results <- do.call(rbind,
                   parallel::mclapply(1:Nrep, function(i) {
                     raw <- compute_both_estimators_SN(M = Ndraws,
                                                       weight = theta)
                     raw$replicate <- i
                     return(raw)
                   }, mc.cores = 8))

results$covers <- is_in(I.quad$value,
                        l = results$lwr,
                        u = results$upr)

library(ggplot2)

ggplot() +
  geom_pointrange(data = results,
                  mapping = aes(
                    x = replicate,
                    y = est,
                    ymin = lwr,
                    ymax = upr
                  )) +
  geom_hline(yintercept = I.quad$value,
             linetype = "dotted") +
  scale_x_continuous("") +
  scale_y_continuous(expression(hat(mu))) +
  facet_grid(proposal ~ ., labeller = label_parsed) +
  theme_bw(base_size = 20)

results

bias_sq <- function(delta) {
  Edelta <- mean(delta)
  return((I.quad$value - Edelta) ^ 2)
}

aggregate(est ~ proposal, var, data = results)
aggregate(est ~ proposal, bias_sq, data = results)

aggregate(var ~ proposal, mean, data = results)
aggregate(covers ~ proposal, mean, data = results)

aggregate(ess ~ proposal, min, data = results)
aggregate(ess ~ proposal, mean, data = results)
aggregate(ess ~ proposal, max, data = results)

subset(results, ess < .5 * Ndraws)

ggplot(results,
       aes(
         x = proposal,
         y = ess,
         colour = proposal,
         fill = proposal
       )) +
  scale_y_continuous("Effective sample size") +
  scale_x_discrete("Proposal",
                   labels = parse(text = unique(results$proposal))) +
  geom_boxplot(alpha = .4) +
  theme_bw(base_size = 16) +
  theme(legend.position = "none")
