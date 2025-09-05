source("evaluation.R")

# --- Plot ----
lambda <- runif(1, 0.5, 10)
nu     <- runif(1, 0.5, 2)
#lambda <- 10
#nu     <- 0.30
log_lambda <- log(lambda)
n_sample <- 4000

# --- Simulate from each method ---
samples_rej_ben <- rcomp_rejection_benson(n_sample, lambda, nu)
samples_rou     <- rcomp_rou(n_sample, lambda, nu)
samples_exa     <- rcomp_exact(n_sample, lambda, nu)
samples_dummy   <- rcomp_dummy(n_sample, lambda, nu)

# --- True pmf (truncate support at reasonable size) ---
K <- max(c(samples_rou, samples_rej_ben, samples_exa, samples_dummy)) + 10
support <- 0:K
log_p <- sapply(support, function(y)
  com_poisson_log_lpmf(y, log_lambda, nu, 1e-12))
p_true <- exp(log_p)
p_true <- p_true / sum(p_true)

df_true <- data.frame(x = support, prob = p_true, type = "True PMF")

# --- Empirical histograms (convert to relative frequency) ---
df_rou     <- data.frame(x = samples_rou, method = "rejection")
df_rej_ben <- data.frame(x = samples_rej_ben, method = "rejection benson")
df_exa     <- data.frame(x = samples_exa, method = "exact")
df_dummy   <- data.frame(x = samples_dummy, method = "dummy")

# Precompute relative frequencies
# RoU sampler
df_rou_freq <- transform(
  as.data.frame(table(df_rou$x)),
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
df_rou_freq$Method <- "Rejection RoU"
df_rej_ben_freq$Method <- "Rejection Benson"
df_exa_freq$Method <- "Exact"
df_dummy_freq$Method <- "Poisson"

# Combine all data
df_plot <- rbind(
  data.frame(x = df_true$x, y = df_true$prob, Method = df_true$Method),
  data.frame(x = df_rou_freq$x, y = df_rou_freq$freq, Method = df_rou_freq$Method),
  data.frame(x = df_rej_ben_freq$x, y = df_rej_ben_freq$freq, Method = df_rej_ben_freq$Method),
  data.frame(x = df_exa_freq$x, y = df_exa_freq$freq, Method = df_exa_freq$Method),
  data.frame(x = df_dummy_freq$x, y = df_dummy_freq$freq, Method = df_dummy_freq$Method)
)

# Plot
library(ggplot2)

ggplot(df_plot, aes(x = x, y = y, color = Method, linetype = Method)) +
  geom_point(size = 2) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c("True PMF" = "black", "Rejection RoU" = "blue", "Rejection Benson" = "orange", "Exact" = "green", "Poisson" = "red")) +
  scale_linetype_manual(values = c("True PMF" = "solid", "Rejection RoU" = "dashed", "Rejection Benson" = "dashed", "Exact" = "dashed", "Poisson" = "dotted")) +
  labs(title = paste0("COM-Poisson λ=", sprintf("%.2f", lambda), ", ν=", sprintf("%.2f", nu)),
       x = "y", y = "Probability / Relative Frequency") +
  theme_minimal()
