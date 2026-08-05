## Project: What makes a green lip? An observational study of age, size, and light effects on the lip colour of farmed greenlip abalone (Haliotis laevigata Donovan)

## Step 1b: Colour correction check (image-level calibration, not per-metric r2)

# Install packages and load libraries
install.packages("here", "ggplot2", "patchwork") 

library(here)
library(ggplot2)
library(patchwork)

# Script 01b: Correction factor structuring, quality/failure summary, and dE spread
### LOAD & INSPECT DATA
df <- read.csv(here("data", "correction_factors.csv"))
nrow(df)
str(df)

##  Calibration quality (good / acceptable / poor) summary 
quality_levels <- c("poor", "acceptable", "good")  # ordered worst to best (matches dE_before_after_by_quality.png in validation/colour_calibration/)
 
# Per-category summary, styled the same way as calibration_summary() above
quality_summary <- function(data, levels = quality_levels) {
  calibrated <- data[data$status == "calibrated", ]
  n_total    <- nrow(calibrated)
  out <- lapply(levels, function(q) {
    n <- sum(calibrated$quality == q, na.rm = TRUE)
    data.frame(
      quality = q,
      n       = n,
      pct     = round(100 * n / n_total, 1),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}
qual_summary <- quality_summary(df)
 
# Print summary
print(qual_summary, row.names = FALSE)
 
## Failed calibration summary (status = "no_checker_found" or "fit_failed") 
fail_statuses <- c("no_checker_found", "fit_failed")
 
failure_summary <- function(data, statuses = fail_statuses) {
  n_total <- nrow(data)
  out <- lapply(statuses, function(s) {
    n <- sum(data$status == s, na.rm = TRUE)
    data.frame(
      status = s,
      n      = n,
      pct    = round(100 * n / n_total, 1),
      stringsAsFactors = FALSE
    )
  })
  base <- do.call(rbind, out)
  n_failed_total <- sum(base$n)
  rbind(base, data.frame(
    status = "all_failed",
    n      = n_failed_total,
    pct    = round(100 * n_failed_total / n_total, 1)
  ))
}
fail_summary <- failure_summary(df)
 
# Print summary
print(fail_summary, row.names = FALSE)
 
### FILTER OUT FAILED CALIBRATION IMAGES 
df_ok <- df[df$status == "calibrated", ]
cat(sprintf("\nRetained %d of %d images (calibration succeeded)\n", nrow(df_ok), nrow(df)))
 
# Order quality as a factor, worst to best, for consistent plotting order
df_ok$quality <- factor(df_ok$quality, levels = quality_levels)
 
### PLOTS: dE spread before calibration, and after calibration by quality 
 
# shared minimal theme: no gridlines, external tick marks, no titles
no_grid_theme <- theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid      = element_blank()
  )
 
# Panel A: dE_before across ALL calibrated images (no quality split) 
mean_before <- mean(df_ok$dE_before, na.rm = TRUE)
sd_before   <- sd(df_ok$dE_before,   na.rm = TRUE)
label_before <- sprintf("%.2f \u00B1 %.2f", mean_before, sd_before)
y_pos_before <- max(df_ok$dE_before, na.rm = TRUE) * 1.08
 
p_before <- ggplot(df_ok, aes(x = "All images", y = dE_before)) +
  geom_boxplot(width = 0.4, fill = "#C44E52", outlier.size = 0.8, outlier.alpha = 0.4) +
  annotate("text", x = 1, y = y_pos_before, label = label_before, size = 3.8) +
  labs(x = NULL, y = expression(Delta * "E")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  no_grid_theme
 
# Panel B: dE_after split by quality (poor / acceptable / good) 
qual_n      <- table(df_ok$quality)
qual_labels <- setNames(paste0(names(qual_n), "\n(n=", as.integer(qual_n), ")"),
                        names(qual_n))
 
# per-quality mean +/- SD of dE_after, positioned just above each group's max
agg_mean <- tapply(df_ok$dE_after, df_ok$quality, mean, na.rm = TRUE)
agg_sd   <- tapply(df_ok$dE_after, df_ok$quality, sd,   na.rm = TRUE)
agg_max  <- tapply(df_ok$dE_after, df_ok$quality, max,  na.rm = TRUE)
 
stats_after <- data.frame(
  quality = factor(names(agg_mean), levels = quality_levels),
  label   = sprintf("%.2f \u00B1 %.2f", agg_mean, agg_sd),
  y_pos   = as.numeric(agg_max) * 1.08
)
 
p_after <- ggplot(df_ok, aes(x = quality, y = dE_after, fill = quality)) +
  geom_boxplot(width = 0.5, outlier.size = 0.8, outlier.alpha = 0.4) +
  geom_text(data = stats_after, aes(x = quality, y = y_pos, label = label),
            inherit.aes = FALSE, size = 3.8) +
  scale_fill_manual(values = c("poor" = "#C44E52", "acceptable" = "#DD8452", "good" = "#55A868")) +
  scale_x_discrete(labels = qual_labels[levels(df_ok$quality)]) +
  labs(x = NULL, y = expression(Delta * "E after calibration")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
  no_grid_theme
 
# Stitch together with patchwork, tagged A / B
combined_plot <- p_before + p_after +
  patchwork::plot_annotation(tag_levels = "A")
 
# Save alongside the other correction-factor outputs
ggsave(filename = here("figures", "dE_before_after_by_quality.png"),plot = combined_plot, width = 10, height = 5, dpi = 300)
 
combined_plot

### END OF SCRIPT ###
