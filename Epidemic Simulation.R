library(MASS)

# ---------------------------------------------------------
# Generate infection-state matrix for a single group
generate_matrix <- function(i, j, mu, alpha, beta, gamma) {
  matrix_result <- matrix(NA, nrow = i, ncol = j)
  # Initial infection state ~ Bernoulli(mu)
  matrix_result[, 1] <- rbinom(i, 1, mu)
  
  for (col in 2:j) {
    I <- sum(matrix_result[, col - 1] == 1)
    for (row in 1:i) {
      if (matrix_result[row, col - 1] == 0) {
        # S → I transition
        p <- 1 - exp(-alpha - beta * I)
      } else {
        # Probability of remaining infected (I → I)
        p <- exp(-gamma)
      }
      matrix_result[row, col] <- rbinom(1, 1, p)
    }
  }
  matrix_result
}

# ---------------------------------------------------------
# Log-likelihood for one group (only α, β, γ)
log_likelihood <- function(params, matrix_result) {
  alpha <- params[1]; beta <- params[2]; gamma <- params[3]
  ll <- 0
  
  for (col in 2:ncol(matrix_result)) {
    I_prev <- sum(matrix_result[, col - 1] == 1)
    p01 <- 1 - exp(-alpha - beta * I_prev)  # S→I
    p00 <- exp(-alpha - beta * I_prev)      # S→S
    p10 <- 1 - exp(-gamma)                  # I→S
    p11 <- exp(-gamma)                      # I→I
    
    for (row in 1:nrow(matrix_result)) {
      s_prev <- matrix_result[row, col - 1]
      s_curr <- matrix_result[row, col]
      if (s_prev == 0 && s_curr == 0) {
        ll <- ll + log(p00)
      } else if (s_prev == 0 && s_curr == 1) {
        ll <- ll + log(p01)
      } else if (s_prev == 1 && s_curr == 0) {
        ll <- ll + log(p10)
      } else if (s_prev == 1 && s_curr == 1) {
        ll <- ll + log(p11)
      }
    }
  }
  ll
}

# ---------------------------------------------------------
# Joint log-prior for α, β, γ  (Gamma(1,1))
joint_prior <- function(params, shape = 1, rate = 1) {
  alpha <- params[1]; beta <- params[2]; gamma <- params[3]
  dgamma(alpha, shape = shape, rate = rate, log = TRUE) +
    dgamma(beta,  shape = shape, rate = rate, log = TRUE) +
    dgamma(gamma, shape = shape, rate = rate, log = TRUE)
}

# ---------------------------------------------------------
# Main algorithm: MH (for α, β, γ) + Gibbs (for μ)
# matrix_data: array [n_individuals, n_times, n_groups]
mh_gibbs_algorithm <- function(matrix_data,
                               initial_params,           
                               proposal_cov,             
                               iterations = 30000,
                               burn_in   = NULL,
                               prior_gamma_shape = 1,
                               prior_gamma_rate  = 1,
                               prior_beta_shape1 = 1,
                               prior_beta_shape2 = 1,
                               adapt = FALSE,            
                               adapt_start = 5000,
                               adapt_rate  = 0.05) {
  if (is.null(burn_in)) burn_in <- round(0.5 * iterations)
  
  n_params <- length(initial_params) # should be 4
  samples  <- matrix(NA, nrow = iterations, ncol = n_params)
  colnames(samples) <- c("alpha","beta","gamma","mu")
  samples[1, ] <- initial_params
  
  n_individuals <- dim(matrix_data)[1]
  n_times       <- dim(matrix_data)[2]
  n_groups      <- dim(matrix_data)[3]
  
  acc <- 0L
  Sigma <- proposal_cov
  
  for (it in 2:iterations) {
    # --- Gibbs update for μ from Beta posterior
    I1 <- sum(matrix_data[, 1, ] == 1)
    n  <- n_individuals * n_groups
    mu_new <- rbeta(1, prior_beta_shape1 + I1, prior_beta_shape2 + (n - I1))
    
    # --- Optional adaptive MH covariance
    if (adapt && it > adapt_start) {
      d <- 3
      s2 <- (2.38^2) / d
      hist_abg <- samples[1:(it-1), 1:3, drop = FALSE]
      if (nrow(hist_abg) > 10) {
        Sigma <- (1 - adapt_rate)^2 * s2 * cov(hist_abg) + 
          (adapt_rate^2) * (0.1^2) * diag(d) / d
      }
    }
    
    # --- MH update for (α, β, γ)
    current_abg <- samples[it - 1, 1:3]
    proposed_abg <- as.numeric(mvrnorm(1, mu = current_abg, Sigma = Sigma))
    
    # Enforce non-negative parameters
    if (any(proposed_abg < 0)) {
      samples[it, ] <- c(current_abg, mu_new)
      next
    }
    
    # Compute log-likelihoods over all groups
    curr_ll <- 0; prop_ll <- 0
    for (g in 1:n_groups) {
      curr_ll <- curr_ll + log_likelihood(current_abg,  matrix_data[ , , g])
      prop_ll <- prop_ll + log_likelihood(proposed_abg, matrix_data[ , , g])
    }
    
    # Add log-priors
    curr_post <- curr_ll + joint_prior(current_abg, prior_gamma_shape, prior_gamma_rate)
    prop_post <- prop_ll + joint_prior(proposed_abg, prior_gamma_shape, prior_gamma_rate)
    
    # MH acceptance ratio
    a_log <- prop_post - curr_post
    if (log(runif(1)) < a_log) {
      current_abg <- proposed_abg
      acc <- acc + 1L
    }
    
    # Save current sample (α, β, γ, μ)
    samples[it, ] <- c(current_abg, mu_new)
  }
  
  list(
    samples = samples,
    accept_rate = acc / (iterations - 1),
    burn_in = burn_in,
    proposal_cov_final = Sigma
  )
}

# =========================================================
# Example: multi-group simulation + MH+Gibbs estimation
# =========================================================

set.seed(123)

# --- True parameter values
alpha_sim <- 0.1
beta_sim  <- 0.1
gamma_sim <- 0.9
mu_sim    <- 0.10

num_individuals <- 10
num_days        <- 100
num_groups      <- 5

# --- Generate 3D simulation data
data_sim <- array(NA, dim = c(num_individuals, num_days, num_groups))
for (g in 1:num_groups) {
  data_sim[, , g] <- generate_matrix(num_individuals, num_days, mu_sim,
                                     alpha_sim, beta_sim, gamma_sim)
}

# --- Priors
prior_gamma_shape <- 1
prior_gamma_rate  <- 1
prior_beta_shape1 <- 1
prior_beta_shape2 <- 1

# --- Initial values and proposal covariance
initial_params <- c(1, 1, 1, 0.5)           # (α, β, γ, μ)
proposal_cov   <- (0.1^2) * diag(3) / 3

# --- Run MH + Gibbs sampler 
iterations <- 30000
burn_in    <- round(0.55 * iterations)

fit <- mh_gibbs_algorithm(
  matrix_data = data_sim,
  initial_params = initial_params,
  proposal_cov = proposal_cov,
  iterations = iterations,
  burn_in = burn_in,
  prior_gamma_shape = prior_gamma_shape,
  prior_gamma_rate  = prior_gamma_rate,
  prior_beta_shape1 = prior_beta_shape1,
  prior_beta_shape2 = prior_beta_shape2,
  adapt = TRUE,
  adapt_start = 5000,
  adapt_rate  = 0.05
)

samples <- fit$samples
cat("Acceptance rate:", round(fit$accept_rate, 3), "\n")

# --- Post burn-in summary
post <- samples[(burn_in + 1):nrow(samples), , drop = FALSE]
theta_hat <- colMeans(post)
print(theta_hat)

# --- Trace plots (showing convergence)
theta_sim <- c(alpha_sim, beta_sim, gamma_sim, mu_sim)
par(mfrow = c(2, 2))
for (k in 1:4) {
  plot(samples[, k], type = "l",
       main = paste0(colnames(samples)[k], " (trace)"),
       xlab = "Iteration", ylab = colnames(samples)[k])
  abline(h = theta_sim[k], col = "red")
}
par(mfrow = c(1, 1))
