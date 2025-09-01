# Unnormalized COM-Poisson mass
log_pi_tilde <- function(x, lambda, nu) {
  x * log(lambda) - nu * lgamma(x + 1)
}

# b_x as in the algorithm
b_x <- function(x, lambda, nu) {
  num <- (x + 1) * sqrt(exp(log_pi_tilde(x + 1, lambda, nu)))
  den <- x * sqrt(exp(log_pi_tilde(x, lambda, nu)))
  return(num / den)
}

# Rejection sampler for n samples
rcomp_rejection <- function(n, lambda, nu) {
  # ---- Setup step ----
  x_star <- floor(lambda^(1 / nu))
  
  while (b_x(x_star, lambda, nu) > 1 && b_x(x_star - 1, lambda, nu) < 1) {
    if (b_x(x_star, lambda, nu) < 1) {
      x_star <- 2 * x_star
    } else {
      x_star <- floor(x_star / 2)
    }
  }
  
  A <- sqrt(exp(log_pi_tilde(floor(lambda^(1 / nu)), lambda, nu)))
  B <- x_star * sqrt(exp(log_pi_tilde(x_star, lambda, nu)))
  
  # ---- Simulation step ----
  out <- numeric(n)
  i <- 1
  while (i <= n) {
    U <- runif(1, 0, A)
    V <- runif(1, 0, B)
    candidate <- V / U
    if (U <= sqrt(exp(log_pi_tilde(floor(candidate), lambda, nu)))) {
      out[i] <- floor(candidate)
      i <- i + 1
    }
  }
  return(out)
}
