source("com_poisson_pmf.R")

# Auxiliar functions
g <- function(y, gamma, nu) {
  if (nu < 1) {
    return(gamma*(1 - gamma)^y)
  } else {
    lresult <- y * log(gamma) - gamma - lgamma(y + 1)
    return(exp(lresult))
  }
}

log_B <- function(mu, nu, p = 1) {
  if (nu < 1) {
    a <- floor(mu / (1 - p)^(1/nu))
    lresult <- (nu * a) * log(mu)
      - a * log(1 - p) - nu * lgamma(a + 1) - log(p)

    return(lresult)
  } else {
    lresult <- (nu - 1) * (floor(mu) * log(mu) - lgamma(floor(mu) + 1))
    
    return(lresult)
  }
}

rcomp_rejection_benson <- function(n, lambda, nu) {
  mu <- lambda^(1/nu)
  out <- numeric(n)
  
  if (nu >= 1) {
    log_b <- log_B(mu, nu)
    
    i <- 1
    while (i <= n) {
      y_sim <- rpois(1, mu)
      lterm <- y_sim * log(mu) - lgamma(y_sim + 1)
      log_alpha <- nu * lterm - log_b - lterm
      u <- runif(1)
      
      if (log(u) < log_alpha) {
        out[i] <- y_sim
        i <- i + 1
      }
    }
  } else {
    p <- (2*nu) / (2*mu*nu + 1 + nu)
    log_b <- log_B(mu, nu, p)
    
    i <- 1
    while (i <= n) {
      y_sim <- rgeom(1, p)
      lterm <- y_sim * log(mu) - lgamma(y_sim + 1)
      log_alpha <- nu * lterm - log_b - lterm
      u <- runif(1, 0, 1)

      if (log(u) < log_alpha) {
        out[i] <- y_sim
        i <- i + 1
      }
    }
  }
  
  return(out)
}

rcomp_rejection_benson(10, 20, 0.1)