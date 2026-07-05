# Project: What makes a green lip? Moderators of lip colour in farmed greenlip abalone (Haliotis laevigata Donovan)

## Step 1: Colour calibration check

# Install packages and load libraries
install.packages("here")
library(here)

# Script 01: Data structuring, factor setup, and inspection

### LOAD & CLEAN DATA
df <- read.csv(here("data", "lip_colour_commercial.csv"))
nrow(df)
str(df)

##  ColorChecker correction quality (r2) per metric

r2_cols <- c("L*" = "L_r2", "a*" = "a_r2", "b*" = "b_r2")  # map each metric to its calibration r2 column 
r2_thresholds <- c(0.95, 0.90, 0.80)  # thresholds 

# Make r2 columns to numeric; treat 0 (and NA) as FAILED corrections 
for (col in r2_cols) {
  df[[col]] <- suppressWarnings(as.numeric(df[[col]]))
}

# Per-metric summary 
calibration_summary <- function(data, cols = r2_cols, thresh = r2_thresholds) {
  out <- lapply(names(cols), function(metric) {
    x      <- data[[cols[[metric]]]]
    failed <- is.na(x) | x == 0            # failed / missing correction
    valid  <- x[!failed]                   # genuine fits only
    base <- data.frame(
      metric    = metric,
      n_total   = length(x),
      n_failed  = sum(failed),
      n_valid   = length(valid),
      mean_r2   = round(mean(valid),   4),
      median_r2 = round(median(valid), 4),
      min_r2    = round(min(valid),    4),
      max_r2    = round(max(valid),    4),
      stringsAsFactors = FALSE
    )
    # % of VALID fits falling below each threshold
    for (t in thresh) {
      base[[paste0("pct_below_", t)]] <- round(100 * mean(valid < t), 1)
    }
    base
  })
  do.call(rbind, out)
}

cal_summary <- calibration_summary(df)

# Print summary
print(cal_summary, row.names = FALSE)

# Breakdown by photographing day (calibration often drifts by session)
calibration_by_day <- function(data, cols = r2_cols, day = "date") {
  do.call(rbind, lapply(names(cols), function(metric) {
    x      <- data[[cols[[metric]]]]
    keep   <- !(is.na(x) | x == 0)
    ag     <- aggregate(x[keep], list(day = data[[day]][keep]),
                        FUN = function(v) c(mean = mean(v), min = min(v),
                                            pct_lt_090 = 100 * mean(v < 0.90)))
    res <- data.frame(metric = metric, day = ag$day,
                      mean_r2    = round(ag$x[, "mean"],       4),
                      min_r2     = round(ag$x[, "min"],        4),
                      pct_r2_lt_090 = round(ag$x[, "pct_lt_090"], 1))
    res
  }))
}

# Print
print(calibration_by_day(df), row.names = FALSE)

# Row-level quality flags 
qc_cut <- 0.90  # flag any < 0.90
for (metric in names(r2_cols)) {
  col <- r2_cols[[metric]]
  ch  <- sub("\\*", "", metric)                      
  df[[paste0("failed_", ch)]] <- is.na(df[[col]]) | df[[col]] == 0
  df[[paste0("poor_",   ch)]] <- !df[[paste0("failed_", ch)]] & df[[col]] < qc_cut
}

cat(sprintf("\nRows flagged poor (valid fit, r2 < %.2f): L*=%d  a*=%d  b*=%d\n",
            qc_cut, sum(df$poor_L), sum(df$poor_a), sum(df$poor_b)))


df$any_bad <- (df$failed_L | df$failed_a | df$failed_b) |
              (df$poor_L   | df$poor_a   | df$poor_b)

# Print date and tank summaries
df$bad_L <- df$failed_L | df$poor_L
df$bad_a <- df$failed_a | df$poor_a
df$bad_b <- df$failed_b | df$poor_b

# Per tank
print(aggregate(cbind(L = bad_L, a = bad_a, b = bad_b) ~ tank, df, sum),
      row.names = FALSE)

# Per date
print(aggregate(cbind(L = bad_L, a = bad_a, b = bad_b) ~ date, df, sum),
      row.names = FALSE)

### END OF SCRIPT ###
