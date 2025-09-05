# Unnormalized COM-Poisson log-mass
log_pi_tilde <- function(x, lambda, nu) {
  x * log(lambda) - nu * lgamma(x + 1)
}

# b_x in log-space
b_x <- function(x, lambda, nu) {
  if (x == 0) {
    return(exp(0.5 * (log_pi_tilde(1, lambda, nu) - log_pi_tilde(0, lambda, nu))))
  }
  log_num <- log(x + 1) + 0.5 * log_pi_tilde(x + 1, lambda, nu)
  log_den <- log(x)     + 0.5 * log_pi_tilde(x, lambda, nu)
  return(exp(log_num - log_den))
}

# Rejection sampler for n samples
rcomp_rou <- function(n, lambda, nu) {
  # ---- Setup step ----
  x_star <- floor(lambda^(1 / nu))
  
  while (b_x(x_star, lambda, nu) > 1 && b_x(x_star - 1, lambda, nu) < 1) {
    if (b_x(x_star, lambda, nu) < 1) {
      x_star <- 2 * x_star
    } else {
      x_star <- floor(x_star / 2)
    }
  }
  
  x0 <- floor(lambda^(1 / nu))
  logA <- 0.5 * log_pi_tilde(x0, lambda, nu)
  A <- exp(logA)
  
  logB <- log(x_star) + 0.5 * log_pi_tilde(x_star, lambda, nu)
  B <- exp(logB)
  
  # ---- Simulation step ----
  out <- numeric(n)
  i <- 1
  while (i <= n) {
    U <- runif(1, 0, A)
    V <- runif(1, 0, B)
    candidate <- round(V / U)
    
    if (log(U) <= 0.5 * log_pi_tilde(candidate, lambda, nu)) {
      out[i] <- candidate
      i <- i + 1
    }
  }
  return(out)
}
