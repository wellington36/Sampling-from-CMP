# --- Setup ---
library(ggplot2)

#set.seed(123)
source("rcomp_rejection.R")   # rejection sampler
source("com_poisson_pmf.R")   # for com_poisson_log_lpmf

# Dummy sampler for test
rcomp_dummy <- function(n, lambda, nu) {
  rpois(n, lambda)  # just Poisson as placeholder
}

# --- Register methods in a list ---
methods <- list(
  rejection = rcomp_rejection,
  dummy     = rcomp_dummy
)

# --- Benchmark wrapper: return multiple metrics ---
benchmark_method <- function(method_fun, n_iter = 20, n_sample = 2000, eps = 1e-12) {
  times  <- numeric(n_iter)
  pvals  <- numeric(n_iter)
  diffs  <- numeric(n_iter)
  
  for (i in seq_len(n_iter)) {
    lambda_i <- runif(1, 1, 10)
    nu_i     <- runif(1, 0.5, 2)
    log_lambda <- log(lambda_i)
    
    # --- timing + sampling ---
    t0 <- Sys.time()
    samples <- method_fun(n = n_sample, lambda = lambda_i, nu = nu_i)
    t1 <- Sys.time()
    times[i] <- as.numeric(difftime(t1, t0, units = "secs"))
    
    # --- sample mean ---
    sample_mean <- mean(samples)
    
    # --- true distribution (truncate at max obs + margin) ---
    K <- max(samples) + 50
    support <- 0:K
    log_p <- sapply(support, function(y)
      com_poisson_log_lpmf(y, log_lambda, nu_i, eps))
    p <- exp(log_p); p <- p / sum(p)
    
    true_mean <- sum(support * p)
    diffs[i] <- abs(sample_mean - true_mean)
    
    # --- chi-square p-value ---
    obs <- table(factor(samples, levels = support))
    pvals[i] <- tryCatch({
      suppressWarnings(chisq.test(obs, p = p)$p.value)
    }, error = function(e) NA)
  }
  
  # Return named metrics
  c(
    Time = mean(times, na.rm = TRUE),
    Chi2Pval = mean(pvals, na.rm = TRUE),
    MeanDiff = mean(diffs, na.rm = TRUE)
  )
}

# --- Run all methods ---
results_list <- lapply(methods, benchmark_method)
results <- do.call(cbind, results_list)
colnames(results) <- names(methods)

# --- Format table for display ---
tab <- rbind(
  c("", colnames(results)),
  c("Time(s)", sprintf("%.4f", results["Time", ])),
  c("Chi2 p-val", sprintf("%.4f", results["Chi2Pval", ])),
  c("|Mean - TrueMean|", sprintf("%.4f", results["MeanDiff", ]))
)

print(tab, quote = FALSE)


# --- Plot ----
lambda <- 1
nu     <- 1.2
log_lambda <- log(lambda)
n_sample <- 5000

# --- Simulate from each method ---
samples_rej   <- rcomp_rejection(n_sample, lambda, nu)
samples_dummy <- rcomp_dummy(n_sample, lambda, nu)

# --- True pmf (truncate support at reasonable size) ---
K <- max(c(samples_rej, samples_dummy)) + 10
support <- 0:K
log_p <- sapply(support, function(y)
  com_poisson_log_lpmf(y, log_lambda, nu, 1e-12))
p_true <- exp(log_p)
p_true <- p_true / sum(p_true)

df_true <- data.frame(x = support, prob = p_true, type = "True PMF")

# --- Empirical histograms (convert to relative frequency) ---
df_rej <- data.frame(x = samples_rej, method = "rejection")
df_dummy <- data.frame(x = samples_dummy, method = "dummy")

# --- Combine into one plot ---
ggplot() +
  geom_col(data = df_true, aes(x = x, y = prob), 
           fill = "grey70", color = "black", width = 0.8, alpha = 0.6) +
  geom_histogram(data = df_rej, aes(x = x, y = ..count../sum(..count..), 
                                    fill = "rejection"), 
                 binwidth = 1, alpha = 0.4, position = "identity") +
  geom_histogram(data = df_dummy, aes(x = x, y = ..count../sum(..count..), 
                                      fill = "dummy"), 
                 binwidth = 1, alpha = 0.4, position = "identity") +
  scale_fill_manual(values = c("rejection" = "blue", "dummy" = "red")) +
  labs(title = paste0("COM-Poisson λ=", lambda, ", ν=", nu),
       x = "y", y = "Probability / Relative Frequency") +
  theme_minimal()
