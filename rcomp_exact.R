source("com_poisson_pmf.R")

# Exact sampler for n samples
rcomp_exact <- function(n, lambda, nu) {
  eps <- 1e-12
  log_lambda <- log(lambda)
  
  # eval blindly until the maximum
  k <- 1
  a <- c(exp(log_k_term(log_lambda, nu, k)))
  
  while (a[k] < exp(log_k_term(log_lambda, nu, k+1)) || k < 2) {
    k <- k + 1
    a <- append(a, exp(log_k_term(log_lambda, nu, k)))
  }
  
  maximum <- k
  
  out <- numeric(n)
  i <- 1
  while (i <= n) {
    z <- runif(1, 0, 1)
    k <- maximum
    j <- 1
    

    F_j <- sum(a[1:j])
    F_k <- sum(a[1:k])
  
    l_k <- F_k
    u_k <- F_k + a[k] * (1 - a[k]/a[k-1])**(-1)
 
    F_j_k_lower <- F_j / u_k
    F_j_k_upper <- F_j / l_k


    while (F_j_k_lower < z) {
      if (F_j_k_lower < z && z < F_j_k_upper) {
        k <- k + 1
      
        if (length(a) < k) {
          a <- append(a, exp(log_k_term(log_lambda, nu, k)))
        }
        F_k <- sum(a[1:k])
      
        l_k <- F_k
        u_k <- F_k + a[k] * (1 - a[k]/a[k-1])**(-1)
      
        F_j_k_lower <- F_j / u_k
        F_j_k_upper <- F_j / l_k
      } else if (z > F_j_k_upper) {
        j <- j + 1
      
        F_j <- sum(a[1:j])
      
        F_j_k_lower <- F_j / u_k
        F_j_k_upper <- F_j / l_k
      }
    }
    out[i] <- j - 1 # -1 because R starts with 1 but we need to start with 0
    i <- i + 1
  }
  return(out)
}

rcomp_exact(1, 2, 0.1)