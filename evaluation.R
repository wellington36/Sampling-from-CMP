# --- Setup ---
library(ggplot2)

#set.seed(121)
source("com_poisson_pmf.R")
source("rcomp_exact.R")       # exact sampler
source("rcomp_rejection.R")   # rejection sampler
source("rcomp_rejection_benson.R") # rejection sampler

# Dummy sampler for test
rcomp_dummy <- function(n, lambda, nu) {
  rpois(n, lambda)  # just Poisson as placeholder
}

# --- Register methods in a list ---
methods <- list(
  rejection        = rcomp_rejection,
  rejection_benson = rcomp_rejection_benson,
  exact            = rcomp_exact,
  dummy            = rcomp_dummy
)

# --- Benchmark wrapper: return multiple metrics ---
benchmark_method <- function(method_fun, n_iter = 20, n_sample = 3000, eps = 1e-10) {
  times  <- numeric(n_iter)
  pvals  <- numeric(n_iter)
  diffs  <- numeric(n_iter)
  
  for (i in seq_len(n_iter)) {
    lambda_i <- runif(1, 0.5, 10)
    nu_i     <- runif(1, 0.5, 2)
    #nu_i <- 1  # Poisson
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
lambda <- 5
nu     <- 2
log_lambda <- log(lambda)
n_sample <- 4000

# --- Simulate from each method ---
samples_rej_ben <- rcomp_rejection_benson(n_sample, lambda, nu)
samples_rej     <- rcomp_rejection(n_sample, lambda, nu)
samples_exa     <- rcomp_exact(n_sample, lambda, nu)
samples_dummy   <- rcomp_dummy(n_sample, lambda, nu)

# --- True pmf (truncate support at reasonable size) ---
K <- max(c(samples_rej, samples_rej_ben, samples_exa, samples_dummy)) + 10
support <- 0:K
log_p <- sapply(support, function(y)
  com_poisson_log_lpmf(y, log_lambda, nu, 1e-12))
p_true <- exp(log_p)
p_true <- p_true / sum(p_true)

df_true <- data.frame(x = support, prob = p_true, type = "True PMF")

# --- Empirical histograms (convert to relative frequency) ---
df_rej     <- data.frame(x = samples_rej, method = "rejection")
df_rej_ben <- data.frame(x = samples_rej_ben, method = "rejection benson")
df_exa     <- data.frame(x = samples_exa, method = "exact")
df_dummy   <- data.frame(x = samples_dummy, method = "dummy")

# Precompute relative frequencies
# Rejection sampler
df_rej_freq <- transform(
  as.data.frame(table(df_rej$x)),
  x = as.numeric(as.character(Var1)),  # <-- convert correctly
  freq = Freq / sum(Freq)
)

# Rejection sampler
df_rej_ben_freq <- transform(
  as.data.frame(table(df_rej_ben$x)),
  x = as.numeric(as.character(Var1)),  # <-- convert correctly
  freq = Freq / sum(Freq)
)

# Exact sampler
df_exa_freq <- transform(
  as.data.frame(table(df_exa$x)),
  x = as.numeric(as.character(Var1)),  # <-- convert correctly
  freq = Freq / sum(Freq)
)

# Dummy sampler
df_dummy_freq <- transform(
  as.data.frame(table(df_dummy$x)),
  x = as.numeric(as.character(Var1)),  # <-- convert correctly
  freq = Freq / sum(Freq)
)


# Add a "Method" column for legend
df_true$Method <- "True PMF"
df_rej_freq$Method <- "Rejection"
df_rej_ben_freq$Method <- "Rejection Benson"
df_exa_freq$Method <- "Exact"
df_dummy_freq$Method <- "Dummy"

# Combine all data
df_plot <- rbind(
  data.frame(x = df_true$x, y = df_true$prob, Method = df_true$Method),
  data.frame(x = df_rej_freq$x, y = df_rej_freq$freq, Method = df_rej_freq$Method),
  data.frame(x = df_rej_ben_freq$x, y = df_rej_ben_freq$freq, Method = df_rej_ben_freq$Method),
  data.frame(x = df_exa_freq$x, y = df_exa_freq$freq, Method = df_exa_freq$Method),
  data.frame(x = df_dummy_freq$x, y = df_dummy_freq$freq, Method = df_dummy_freq$Method)
)

# Plot
library(ggplot2)

ggplot(df_plot, aes(x = x, y = y, color = Method, linetype = Method)) +
  geom_point(size = 2) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("True PMF" = "black", "Rejection" = "blue", "Rejection Benson" = "orange", "Exact" = "green", "Dummy" = "red")) +
  scale_linetype_manual(values = c("True PMF" = "solid", "Rejection" = "dashed", "Rejection Benson" = "dashed", "Exact" = "dashed", "Dummy" = "dotted")) +
  labs(title = paste0("COM-Poisson λ=", lambda, ", ν=", nu),
       x = "y", y = "Probability / Relative Frequency") +
  theme_minimal()
