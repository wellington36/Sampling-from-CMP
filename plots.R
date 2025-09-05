source("evaluation.R")

# --- Choose methods to use ---
methods_to_use <- c("rou", "rej_ben", "exa", "dummy")  
# options: "rou", "rej_ben", "exa", "dummy"

# --- Parameters ---
lambda <- runif(1, 0.5, 10)
nu     <- runif(1, 0.5, 2)
log_lambda <- log(lambda)
n_sample <- 4000

# --- Storage ---
samples_list <- list()

# --- Simulate conditionally ---
if ("rou" %in% methods_to_use) {
  samples_list$rou <- rcomp_rou(n_sample, lambda, nu)
}
if ("rej_ben" %in% methods_to_use) {
  samples_list$rej_ben <- rcomp_rejection_benson(n_sample, lambda, nu)
}
if ("exa" %in% methods_to_use) {
  samples_list$exa <- rcomp_exact(n_sample, lambda, nu)
}
if ("dummy" %in% methods_to_use) {
  samples_list$dummy <- rcomp_dummy(n_sample, lambda, nu)
}

# --- True pmf (truncate support at reasonable size) ---
K <- max(unlist(samples_list)) + 10
support <- 0:K
log_p <- sapply(support, function(y)
  com_poisson_log_lpmf(y, log_lambda, nu, 1e-12))
p_true <- exp(log_p)
p_true <- p_true / sum(p_true)

df_true <- data.frame(
  x = support,
  y = p_true,        # use 'y' instead of 'prob'
  Method = "True PMF"
)


# --- Empirical histograms ---
make_freq_df <- function(samples, name) {
  tbl <- as.data.frame(table(x = samples))  # explicitly name the column 'x'
  data.frame(
    x = as.numeric(as.character(tbl$x)),
    y = tbl$Freq / sum(tbl$Freq),
    Method = name
  )
}

df_plot <- list(df_true)

if ("rou" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$rou, "Rejection RoU")))
}
if ("rej_ben" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$rej_ben, "Rejection Benson")))
}
if ("exa" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$exa, "Exact")))
}
if ("dummy" %in% methods_to_use) {
  df_plot <- c(df_plot, list(make_freq_df(samples_list$dummy, "Poisson")))
}

df_plot <- do.call(rbind, df_plot)

# --- Plot ---
library(ggplot2)

ggplot(df_plot, aes(x = x, y = y, color = Method, linetype = Method)) +
  geom_point(size = 2) +
  geom_line(size = 0.5) +
  scale_color_manual(values = c(
    "True PMF" = "black",
    "Rejection RoU" = "blue",
    "Rejection Benson" = "orange",
    "Exact" = "green",
    "Poisson" = "red"
  )) +
  scale_linetype_manual(values = c(
    "True PMF" = "solid",
    "Rejection RoU" = "dashed",
    "Rejection Benson" = "dashed",
    "Exact" = "dashed",
    "Poisson" = "dotted"
  )) +
  labs(title = paste0("COM-Poisson λ=", sprintf("%.2f", lambda), 
                      ", ν=", sprintf("%.2f", nu)),
       x = "y", y = "Probability / Relative Frequency") +
  theme_minimal()
