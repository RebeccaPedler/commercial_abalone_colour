# Project: What makes a green lip? Moderators of lip colour in farmed greenlip abalone (Haliotis laevigata Donovan)

## Step 2: Running linear mixed effect models

# Install packages and load libraries

install.packages(c("tidyverse", "ggplot2", "patchwork", "lme4", "lmerTest", "emmeans", "performance", "nlme", "ggcorrplot", "farver"))

library(tidyverse)
library(ggplot2)
library(patchwork)
library(lme4)        
library(lmerTest)    
library(emmeans)     
library(performance) 
library(car)  
library(tidyverse)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)    
library(nlme)
library(ggcorrplot)
library(farver)

# Create directory to save all figures 
fig_dir <- here::here("figures")

# Script 01: Data structuring, factor setup, and inspection

### LOAD & CLEAN DATA
df_raw <- read.csv(here("data", "lip_colour_commercial.csv"))
nrow(df_raw)
str(df_raw)

# Rename columns
names(df_raw)[names(df_raw) == "lightness"] <- "raw_lightness"
names(df_raw)[names(df_raw) == "a"] <- "raw_A"
names(df_raw)[names(df_raw) == "b"] <- "raw_B"

names(df_raw)[names(df_raw) == "lightness_corrected"] <- "lightness"
names(df_raw)[names(df_raw) == "a_corrected"] <- "a"
names(df_raw)[names(df_raw) == "b_corrected"] <- "b"

df_raw <- df_raw |>
  mutate(across(where(is.character), str_trim)) |>
  filter(!is.na(farm)) |>
  filter(lightness != "#N/A") |>
  filter(raw_lightness != 0) |>
  filter(lightness != 0) 

# Correct data types
df <- df_raw |>
  mutate(
    date = factor(date),

    # Unordered factors
    section_coverage = factor(section_coverage),
    farm             = factor(farm),
    section          = factor(section),
    tank             = factor(tank),
    age              = factor(age, levels = c(19, 31)),
    diet             = factor(diet),

    # Numeric
    lightness        = as.numeric(lightness),
    a                = as.numeric(a),
    b                = as.numeric(b),
    chroma           = as.numeric(chroma),
    length           = as.numeric(length_mm),
    width            = as.numeric(width_mm),
    area             = as.numeric(area_mm2),

    abalone_id = paste(tank, abalone_number, sep = "_")
  )

### CLEAN DATA 

# Flag duplicate image IDs
dupes <- df_raw |> filter(duplicated(image_ID) | duplicated(image_ID, fromLast = TRUE))
if (nrow(dupes) > 0) {
  cat("WARNING: Duplicate image IDs found:\n")
  print(dupes[, c("image_ID", "tank", "abalone_number", "length_mm")], row.names = FALSE)
} else {
  cat("No duplicate image IDs.\n")
}

# Check ranges for CIELAB values and length
cat(sprintf("  L*: min=%.2f, max=%.2f, median=%.2f\n",
            min(df$lightness, na.rm=TRUE), max(df$lightness, na.rm=TRUE), median(df$lightness, na.rm=TRUE))) # should be (0 - 100)
cat(sprintf("  a*: min=%.2f, max=%.2f, median=%.2f\n",
            min(df$a, na.rm=TRUE), max(df$a, na.rm=TRUE), median(df$a, na.rm=TRUE))) 
cat(sprintf("  b*: min=%.2f, max=%.2f, median=%.2f\n",
            min(df$b, na.rm=TRUE), max(df$b, na.rm=TRUE), median(df$b, na.rm=TRUE)))

# Summary data for abalone morphometrics
df_no_length_zero <- df_raw |> filter(length_mm != 0)

cat(sprintf("  length_mm*: min=%.2f, max=%.2f, median=%.2f\n",
            min(df_no_length_zero$length_mm, na.rm=TRUE), max(df_no_length_zero$length_mm, na.rm=TRUE), median(df_no_length_zero$length_mm, na.rm=TRUE)))

cat(sprintf("  width_mm*: min=%.2f, max=%.2f, median=%.2f\n",
            min(df_no_length_zero$width_mm, na.rm=TRUE), max(df_no_length_zero$width_mm, na.rm=TRUE), median(df_no_length_zero$width_mm, na.rm=TRUE)))

cat(sprintf("  area_mm2*: min=%.2f, max=%.2f, median=%.2f\n",
            min(df_no_length_zero$area_mm2, na.rm=TRUE), max(df_no_length_zero$area_mm2, na.rm=TRUE), median(df_no_length_zero$area_mm2, na.rm=TRUE)))

# Flag any L* outside valid range (if nothing is returned, there are no outliers)
invalid_L <- df |> filter(lightness < 0 | lightness > 100)
if (nrow(invalid_L) > 0) {
  cat("WARNING: Invalid L* values:\n")
  print(invalid_L[, c("image_ID", "tank", "lightness")], row.names = FALSE)
}

# Flag within-tank outliers (>3 SD ) 
flag_outliers <- function(data, var) {
  data |>
    group_by(tank) |>
    mutate(
      tank_mean = mean(.data[[var]], na.rm = TRUE),
      tank_sd   = sd(.data[[var]],   na.rm = TRUE),
      z_within  = (.data[[var]] - tank_mean) / tank_sd,
      outlier   = abs(z_within) > 3
    ) |>
    ungroup() |>
    filter(outlier) |>
    dplyr::select(image_ID, tank, diet, age, all_of(var), z_within)
}

for (var in c("lightness", "a", "b")) {
  out <- flag_outliers(df, var)
  cat(sprintf("\n  Within-tank outliers (|z|>3) for %s: %d rows\n", var, nrow(out)))
  if (nrow(out) > 0) print(as.data.frame(out), row.names = FALSE)
}

### OBSERVE LEVELS AND OBSERVATIONS FOR EACH FACTOR

# Total number of observations
nrow(df)

factor_vars <- c("farm", "section", "age", "diet",
                 "section_coverage")

for (v in factor_vars) {
  cat(sprintf("\n--- %s ---\n", v))
  df |>
    count(.data[[v]], name = "n") |>
    mutate(pct = round(n / sum(n) * 100, 1)) |>
    as.data.frame() |>
    print(row.names = FALSE)
}

# Tank-level summaries
# Tanks per section_coverage
df |>
  distinct(tank, section_coverage) |>
  count(section_coverage) |>
  as.data.frame() |>
  print(row.names = FALSE)

# Abalone per tank
df |>
  count(tank, name = "n_abalone") |>
  summarise(mean = round(mean(n_abalone), 1),
            sd   = round(sd(n_abalone),   1),
            min  = min(n_abalone),
            max  = max(n_abalone)) |>
  as.data.frame() |>
  print(row.names = FALSE)

### GENERAL DATA OBSERVATIONS BY PLOTTING

# Create colour palettes for plotting
sec_cols <- c("single" = "#8DB4C8", "double" = "#2E6E9E") 
age_cols  <- c("19"     = "#E07B39", "31"     = "#5B8FA8")   
 
## Overall distributions and histograms — L*, a*, b*
hist_plot <- function(var, xlab, fill_col) {
  ggplot(df, aes(x = .data[[var]])) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 40, fill = fill_col,
                   colour = "white", linewidth = 0.25, alpha = 0.85) +
    geom_density(colour = "#1B2A1C", linewidth = 0.6) +
    geom_rug(colour = fill_col, alpha = 0.3, linewidth = 0.3) +
    labs(x = xlab, y = "Density") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          axis.title.y     = element_text(size = 9, colour = "grey40"))
}

p_L <- hist_plot("lightness", "Lightness",   "#5B8FA8")
p_a <- hist_plot("a", "a* (green\u2013red)",   "#7A9E5A")
p_b <- hist_plot("b", "b* (blue\u2013yellow)", "#C4A24A")
p_chroma <- hist_plot("chroma", "chroma", "#E67E22")

# Patch plots
p_hist_all <- (p_L + ggtitle("A)") +
                 theme(plot.title = element_text(face = "bold"))) |
              (p_a + ggtitle("B)") +
                 theme(plot.title = element_text(face = "bold"))) |
              (p_b + ggtitle("C)") +
                 theme(plot.title = element_text(face = "bold")))

print(p_hist_all)

ggsave(file.path(fig_dir, "p_hist_all.png"), plot = p_hist_all, dpi = 300, width = 12, height = 8, units = "in")

# Print chroma 
print(p_chroma)

ggsave(file.path(fig_dir, "p_chroma.png"), plot = p_chroma, dpi = 300, width = 8, height = 8, units = "in")

## Inspect section_coverage (single or double shade cloth)

lab_long <- function(data) {
  data |>
    pivot_longer(cols = c(lightness, a, b),
                 names_to = "metric", values_to = "value") |>
    mutate(metric = factor(metric,
                           levels = c("lightness", "a", "b"),
                           labels = c("Lightness", "a*", "b*")))
}

# Histograms faceted by section_coverage
p_sec_hist <- df |>
  lab_long() |>
  ggplot(aes(x = value, fill = section_coverage)) +
  geom_histogram(bins = 30, colour = "white", linewidth = 0.2,
                 alpha = 0.8, position = "identity") +
  facet_grid(section_coverage ~ metric, scales = "free_x") +
  scale_fill_manual(values = sec_cols, guide = "none") +
  labs(x = "Value", y = "Count") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line = element_line(colour = "black"),  
    strip.text = element_text(size = 10),
    plot.title = element_text(size = 13, face = "bold"),
    plot.subtitle = element_text(size = 10, colour = "grey40")
  )

# Print plot
print(p_sec_hist)
ggsave(file.path(fig_dir, "p_sec_hist.png"), plot = p_sec_hist, dpi = 300, width = 12, height = 8, units = "in")

# Boxplots by section_coverage (no age)
p_sec_box <- df |>
  lab_long() |>
  ggplot(aes(x = section_coverage, y = value, fill = section_coverage)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5,
               colour = "grey30", linewidth = 0.4, alpha = 0.8) +
  geom_jitter(width = 0.15, alpha = 0.07, size = 0.6, colour = "grey20") +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = sec_cols, guide = "none") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 11),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"))

# Print plot
print(p_sec_box)
ggsave(file.path(fig_dir, "p_sec_box.png"), plot = p_sec_box, dpi = 300, width = 12, height = 8, units = "in")

## Inspect section_coverage (single or double shade cloth) X age (19 and 31)

# Boxplots: section_coverage × age (dodged)
p_sec_age_box <- df |>
  lab_long() |>
  ggplot(aes(x = section_coverage, y = value, fill = age)) +
  geom_boxplot(outlier.size = 1, outlier.alpha = 0.5,
               colour = "grey30", linewidth = 0.4, alpha = 0.8,
               position = position_dodge(0.8)) +
  facet_wrap(~metric, scales = "free_y") +
  scale_fill_manual(values = age_cols, name = "Age (months)") +
  labs(x = "Section coverage", y = "Value") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 11),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"))

print(p_sec_age_box)
ggsave(file.path(fig_dir, "p_sec_age_box.png"), plot = p_sec_age_box, dpi = 300, width = 12, height = 8, units = "in")

# Histograms: section_coverage × age
p_sec_age_hist <- df |>
  lab_long() |>
  ggplot(aes(x = value, fill = age)) +
  geom_histogram(bins = 25, colour = "white", linewidth = 0.2,
                 alpha = 0.75, position = "identity") +
  facet_grid(section_coverage ~ metric, scales = "free_x") +
  scale_fill_manual(values = age_cols, name = "Age (months)") +
  labs(x = "Value", y = "Count")+
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 11),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 13, face = "bold"),
        plot.subtitle = element_text(size = 10, colour = "grey40"))

# Print plot
print(p_sec_age_hist)
ggsave(file.path(fig_dir, "p_sec_age_hist.png"), plot = p_sec_age_hist, dpi = 300, width = 12, height = 8, units = "in")

# Create numeric summaries
# By section coverage
df |>
  group_by(section_coverage) |>
  summarise(n           = n(),
            L_mean      = round(mean(lightness), 2), L_sd      = round(sd(lightness), 2),
            a_mean      = round(mean(a),         2), a_sd      = round(sd(a),         2),
            b_mean      = round(mean(b),         2), b_sd      = round(sd(b),         2),
            chroma_mean = round(mean(chroma),    2), chroma_sd = round(sd(chroma),    2),
            .groups = "drop") |>
  as.data.frame() |>
  print(row.names = FALSE)

# By age
df |>
  group_by(age) |>
  summarise(n           = n(),
            L_mean      = round(mean(lightness), 2), L_sd      = round(sd(lightness), 2),
            a_mean      = round(mean(a),         2), a_sd      = round(sd(a),         2),
            b_mean      = round(mean(b),         2), b_sd      = round(sd(b),         2),
            chroma_mean = round(mean(chroma),    2), chroma_sd = round(sd(chroma),    2),
            .groups = "drop") |>
  as.data.frame() |>
  print(row.names = FALSE)

# By age and section_coverage
df |>
  group_by(section_coverage, age) |>
  summarise(n           = n(),
            L_mean      = round(mean(lightness), 2), L_sd      = round(sd(lightness), 2),
            a_mean      = round(mean(a),         2), a_sd      = round(sd(a),         2),
            b_mean      = round(mean(b),         2), b_sd      = round(sd(b),         2),
            chroma_mean = round(mean(chroma),    2), chroma_sd = round(sd(chroma),    2),
            .groups = "drop") |>
  as.data.frame() |>
  print(row.names = FALSE)

### Script 02: Collinearity assessment — section_coverage, age, diet

# Create separate dataframe to work on tank level
tank_df <- df |>
  distinct(tank, section_coverage, age, diet, farm, section)

## CROSSTAB COUNTS (tank level)

# section_coverage x age
table(tank_df$section_coverage, tank_df$age) |>
  addmargins() |>
  print()

# section_coverage
table(tank_df$section_coverage, tank_df$diet) |>
  addmargins() |>
  print()

# age x diet
table(tank_df$age, tank_df$diet) |>
  addmargins() |>
  print()

# Three-way: section_coverage x age x diet (n tanks)
tank_df |>
  count(section_coverage, age, diet) |>
  arrange(section_coverage, age, diet) |>
  as.data.frame() |>
  print(row.names = FALSE)

## VISUALISE 

p_tile <- tank_df |>
  count(section_coverage, age, diet) |>
  ggplot(aes(x = age, y = diet, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = n), size = 3.5, fontface = "bold",
            colour = "white") +
  facet_wrap(~section_coverage, labeller = label_both) +
  scale_fill_gradient(low = "#BDD7EE", high = "#1F618D",
                      name = "Tanks (n)") +
  labs(x = "Age", y = "Diet") +
  theme_minimal(base_size = 11) +
  theme(panel.grid    = element_blank(),
        strip.text    = element_text(size = 11, face = "bold"),
        plot.title    = element_text(size = 13, face = "bold"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.line = element_line(colour = "black"))

# Print crossing plot
print(p_tile)

# Save plots
ggsave(file.path(fig_dir, "p_tile.png"), plot = p_tile, width = 8, height = 8, dpi = 300)

## PAIRWISE CROSSTABS WITH CRAMER'S V

cramer_v <- function(x, y) {
  tbl <- table(x, y)
  chi <- suppressWarnings(chisq.test(tbl))
  n   <- sum(tbl)
  k   <- min(nrow(tbl), ncol(tbl))
  v   <- sqrt(chi$statistic / (n * (k - 1)))
  p   <- chi$p.value
  list(V = round(as.numeric(v), 3), p = round(p, 4),
       chi2 = round(chi$statistic, 2), df = chi$parameter)
}

# Create pairs for moderators
pairs <- list(
  c("section_coverage", "age"),
  c("section_coverage", "diet"),
  c("age",              "diet")
)

for (pair in pairs) {
  res <- cramer_v(tank_df[[pair[1]]], tank_df[[pair[2]]])
  cat(sprintf("%-20s x %-20s  V = %.3f  p = %.4f\n",
              pair[1], pair[2], res$V, res$p))
}

## VARIANCE INFLATION — manual check via linear probability models

tank_num <- tank_df |>
  mutate(
    sec_double = as.integer(section_coverage == "double"),
    age_19      = as.integer(age == "19"),
    # diet reference = 1
    diet_2 = as.integer(diet == 2),
    diet_3 = as.integer(diet == 3),
    diet_4 = as.integer(diet == 4),
    diet_5 = as.integer(diet == 5)
  )

calc_vif <- function(outcome, predictors, data) {
  f   <- as.formula(paste(outcome, "~", paste(predictors, collapse = " + ")))
  mod <- lm(f, data = data)
  r2  <- summary(mod)$r.squared
  vif <- 1 / (1 - r2)
  cat(sprintf("  VIF for %-14s = %.2f  (R² = %.3f)\n", outcome, vif, r2
  ))
}

calc_vif("sec_double", c("age_19","diet_2","diet_3","diet_4","diet_5"), tank_num)
calc_vif("age_19",c("sec_double","diet_2","diet_3","diet_4","diet_5"),tank_num)
calc_vif("diet_2", c("sec_double","age_19"), tank_num)
calc_vif("diet_3", c("sec_double","age_19"), tank_num)
calc_vif("diet_4", c("sec_double","age_19"), tank_num)
calc_vif("diet_5", c("sec_double","age_19"), tank_num)

# Diets not evenly represented accross age.

## CORRELATION OF CIELAB METRICS

# Prepare data
colour_df <- df |>
  dplyr::select(L = lightness, a, b, C = chroma) |>
  tidyr::drop_na()
 
# Pearson correlation + p-values
cor_r <- cor(colour_df, method = "pearson")
cor_p <- cor_pmat(colour_df)
 
# Print results
print(round(cor_r, 3))
print(format(cor_p, scientific = TRUE))
 
# Plot
CIELAB_correlation_p <- ggcorrplot(cor_r,
                p.mat      = cor_p,
                method     = "square",
                type       = "lower",
                lab        = TRUE,
                lab_size   = 5,
                sig.level  = 0.05,
                insig      = "blank",
                colors     = c("#C0392B", "white", "#1F618D")) +        
        theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
 
print(CIELAB_correlation_p)
ggsave(file.path(fig_dir, "CIELAB_correlation_p.png"), CIELAB_correlation_p, width = 10, height = 10, dpi = 300)

### Script 03: Linear Mixed Models — b*, L*, a* ~ section_coverage + age

## PREPARE MODELLING DATASET

df_mod <- df |>
  filter(!is.na(section_coverage), !is.na(age),
         !is.na(b), !is.na(a), !is.na(lightness)) |>
  mutate(
    farm             = factor(farm, ordered = FALSE),
    section_coverage = factor(section_coverage, ordered = FALSE),
    diet             = factor(diet, ordered = FALSE),
    age              = factor(age, ordered = FALSE),
    section_coverage = relevel(section_coverage, ref = "single"),
    age              = relevel(age, ref = "19"),
    diet             = relevel(diet, ref = "1")
  )

## NULL MODELS & ICC 
# Fit intercept-only models (REML) to estimate ICC before any fixed effects

null_b <- lmer(b ~ 1 + (1 | tank), data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))
null_L <- lmer(lightness ~ 1 + (1 | tank), data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))
null_a <- lmer(a ~ 1 + (1 | tank), data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))
null_chroma <- lmer(chroma ~ 1 + (1 | tank), data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))

icc_table <- data.frame(
  response = c("b*", "L*", "a*", "chroma"),
  ICC      = round(c(performance::icc(null_b)$ICC_adjusted,
                     performance::icc(null_L)$ICC_adjusted,
                     performance::icc(null_a)$ICC_adjusted,
                     performance::icc(null_chroma)$ICC_adjusted), 4
))

print(icc_table, row.names = FALSE)

## FORWARD SELECTION: FIT ALL MODEL COMBINATIONS AND COMPARE

# Order: null → + section_coverage → + age (in order of hypothetic importance)

fit_steps <- function(response, data = df_mod) {
  f0  <- as.formula(paste(response, "~ 1                                   + (1 | tank)"))
  f1  <- as.formula(paste(response, "~ section_coverage                    + (1 | tank)"))
  f2  <- as.formula(paste(response, "~ age                                 + (1 | tank)"))
  f3  <- as.formula(paste(response, "~ section_coverage + age              + (1 | tank)"))

  list(
    m0 = lmer(f0, data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa")),
    m1 = lmer(f1, data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa")),
    m2 = lmer(f2, data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa")),
    m3 = lmer(f3, data = data, REML = FALSE, control = lmerControl(optimizer = "bobyqa"))
  )
}

models_b <- fit_steps("b")
models_L <- fit_steps("lightness")
models_a <- fit_steps("a")
models_chroma <- fit_steps("chroma")

# AIC / BIC COMPARISON TABLE 

model_labels <- c(
  "m0: null (intercept only)",
  "m1: + section_coverage",
  "m2: + age",
  "m3: + section_coverage + age"
)

aic_table <- function(mods, response) {
  aic_vals <- sapply(mods, AIC)
  bic_vals <- sapply(mods, BIC)
  loglik   <- sapply(mods, logLik)
  df_used  <- sapply(mods, function(m) attr(logLik(m), "df"))

  data.frame(
    response  = response,
    model     = model_labels,
    df        = df_used,
    logLik    = round(loglik,   2),
    AIC       = round(aic_vals, 2),
    dAIC      = round(aic_vals - min(aic_vals), 2),
    BIC       = round(bic_vals, 2),
    dBIC      = round(bic_vals - min(bic_vals), 2)
  )
}

aic_all <- bind_rows(
  aic_table(models_L, "lightness"),
  aic_table(models_a, "a*"),
  aic_table(models_b, "b*"),
  aic_table(models_chroma, "chroma")
)

# Print per response 
print(aic_all)

# Refit final models 

## REFIT SELECTED MODEL WITH REML

final_formula_L      <- lightness ~ section_coverage + age + (1 | tank)
final_formula_a      <- a         ~ section_coverage + age + (1 | tank)
final_formula_b      <- b         ~ section_coverage + age + (1 | tank)
final_formula_chroma <- chroma    ~ section_coverage + age + (1 | tank)

full_L <- lmer(final_formula_L, data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))
full_a <- lmer(final_formula_a, data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))
full_b <- lmer(final_formula_b, data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))
full_chroma <- lmer(final_formula_chroma, data = df_mod, REML = TRUE,
               control = lmerControl(optimizer = "bobyqa"))

summary(full_L)
summary(full_a)
summary(full_b)
summary(full_chroma)

### Script 04: LMM diagnostics — Lightness, a*, b*, chroma

# MODEL LIST 
model_list <- list(
  list(mod = full_L, label = "L*",     response = "lightness"),
  list(mod = full_a, label = "a*",     response = "a"),
  list(mod = full_b, label = "b*",     response = "b"),
  list(mod = full_chroma, label = "Chroma", response = "chroma")
)

# DIAGNOSTIC FUNCTION 

run_diagnostics <- function(mod, label, resp, data = df_mod) {

  fitted_vals <- fitted(mod)
  resid_raw   <- residuals(mod, type = "pearson")
  resid_sc    <- scale(resid_raw)[, 1]
  blups       <- ranef(mod)$tank[, 1]

  # Shapiro-Wilk on residuals
  set.seed(42)
  sw_r <- shapiro.test(sample(resid_sc, min(length(resid_sc), 4999)))

  # Shapiro-Wilk on BLUPs
  sw_b <- shapiro.test(blups)

  # Levene-type test
  df_test <- data |> mutate(resid_abs = abs(resid_raw))

  for (grp in c("section_coverage", "age")) {
    fit <- lm(as.formula(paste("resid_abs ~", grp)), data = df_test)
    fst <- summary(fit)$fstatistic
    p   <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)
    sig <- case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                     p < 0.05  ~ "*",   p < 0.1  ~ ".",
                     TRUE ~ "ns")
  }

  # Kruskal-Wallis: BLUP independence
  blup_df <- data.frame(blup = blups,
                        tank = rownames(ranef(mod)$tank)) |>
    left_join(
      data |> distinct(tank, section_coverage, age) |>
              mutate(tank = as.character(tank)),
      by = "tank")

  for (grp in c("section_coverage", "age")) {
    kw  <- kruskal.test(as.formula(paste("blup ~", grp)), data = blup_df)
    sig <- case_when(kw$p.value < 0.001 ~ "***",
                     kw$p.value < 0.05  ~ "*",
                     TRUE               ~ "ns")
  }

  # Residual ICC — FIX: definition was broken by an inline comment
  resid_df    <- data.frame(tank = data$tank, resid = resid_raw)
  overall_var <- var(resid_df$resid)
  mean_within <- resid_df |>
                   group_by(tank) |>
                   summarise(v = var(resid), .groups = "drop") |>
                   pull(v) |> mean()
  resid_icc   <- max(0, 1 - mean_within / overall_var)

# Shapiro-Wilk on residuals
  set.seed(42)
  sw_r <- shapiro.test(sample(resid_sc, min(length(resid_sc), 4999)))
  cat(sprintf("  SW residuals:  W = %.4f, p = %.4f\n", sw_r$statistic, sw_r$p.value))

  # Shapiro-Wilk on BLUPs
  sw_b <- shapiro.test(blups)
  cat(sprintf("  SW BLUPs:      W = %.4f, p = %.4f\n", sw_b$statistic, sw_b$p.value))

  # Levene-type test
  df_test <- data |> mutate(resid_abs = abs(resid_raw))
  for (grp in c("section_coverage", "age")) {
    fit <- lm(as.formula(paste("resid_abs ~", grp)), data = df_test)
    fst <- summary(fit)$fstatistic
    p   <- pf(fst[1], fst[2], fst[3], lower.tail = FALSE)
    sig <- case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                     p < 0.05  ~ "*",   p < 0.1  ~ ".", TRUE ~ "ns")
    cat(sprintf("  Levene %-20s F = %.3f, p = %.4f %s\n", grp, fst[1], p, sig))
  }

  # Kruskal-Wallis
  for (grp in c("section_coverage", "age")) {
    kw  <- kruskal.test(as.formula(paste("blup ~", grp)), data = blup_df)
    sig <- case_when(kw$p.value < 0.001 ~ "***", kw$p.value < 0.05 ~ "*", TRUE ~ "ns")
    cat(sprintf("  KW BLUP %-20s chi2 = %.3f, p = %.4f %s\n", grp, kw$statistic, kw$p.value, sig))
  }

  # Residual ICC
  cat(sprintf("  Residual ICC:  %.4f\n", resid_icc))

  # PLOTS
  df_p <- data.frame(fitted = fitted_vals,
                     resid  = resid_raw,
                     sq_abs = sqrt(abs(resid_sc)))

  p1 <- ggplot(df_p, aes(x = fitted, y = resid)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.15, size = 0.7, colour = "#2E6E9E") +
    geom_smooth(method = "loess", se = TRUE, span = 0.6,
                colour = "#E07B39", fill = "#E07B39",
                alpha = 0.12, linewidth = 0.8) +
    labs(x = "Fitted", y = "Pearson residuals",
         title = paste0("Residuals vs fitted — ", label),
         subtitle = "Linearity + homoscedasticity: loess should be flat, spread even") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey40"))

  p2 <- ggplot(data.frame(r = resid_sc), aes(sample = r)) +
    stat_qq(alpha = 0.2, size = 0.7, colour = "#2E6E9E") +
    stat_qq_line(colour = "#E07B39", linewidth = 0.8) +
    labs(x = "Theoretical quantiles", y = "Scaled residuals",
         title = paste0("QQ plot — ", label),
         subtitle = "Normality: points should follow the line") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey40"))

  p3 <- ggplot(df_p, aes(x = fitted, y = sq_abs)) +
    geom_point(alpha = 0.15, size = 0.7, colour = "#2E6E9E") +
    geom_smooth(method = "loess", se = TRUE, span = 0.6,
                colour = "#E07B39", fill = "#E07B39",
                alpha = 0.12, linewidth = 0.8) +
    labs(x = "Fitted", y = expression(sqrt("|Scaled residuals|")),
         title = paste0("Scale-location — ", label),
         subtitle = "Homoscedasticity: loess should be flat") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey40"))

  p4 <- ggplot(data.frame(b = blups), aes(sample = b)) +
    stat_qq(colour = "#5B8FA8", size = 2, alpha = 0.85) +
    stat_qq_line(colour = "#E07B39", linewidth = 0.8) +
    labs(x = "Theoretical quantiles", y = "BLUPs",
         title = paste0("BLUP QQ — ", label),
         subtitle = paste0("n = ", length(blups), " tanks")) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey40"))

  panel <- (p1 | p2 | p3 | p4) +
    plot_annotation(
      title = paste0("Diagnostic plots — ", label),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )

  print(panel)
  ggsave(file.path(fig_dir, sprintf("plot_diagnostics_%s.png", resp)),
         panel, width = 16, height = 4, dpi = 300)
  cat(sprintf("\n  Saved: plot_diagnostics_%s.png\n", resp))

} 

# RUN FOR ALL MODELS 
run_diagnostics(mod = full_L, label = "L*", resp = "L")
run_diagnostics(mod = full_a, label = "a*", resp = "a")
run_diagnostics(mod = full_b, label = "b*", resp = "b")
run_diagnostics(mod = full_chroma, label = "chroma*", resp = "chroma")

# Script 05: Weighted LME models (varIdent) 

## FIT WEIGHTED MODELS USING VARiIDENT
#    varIdent(form = ~1|section_coverage) - most heterogeneous

ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

weighted_L <- lme(
  lightness      ~ section_coverage + age,
  random         = ~ 1 | tank,
  weights        = varIdent(form = ~ 1 | section_coverage),
  data           = df_mod,
  method         = "REML",
  control        = ctrl
)

weighted_a <- lme(
  a              ~ section_coverage + age,
  random         = ~ 1 | tank,
  weights        = varIdent(form = ~ 1 | section_coverage),
  data           = df_mod,
  method         = "REML",
  control        = ctrl
)

weighted_b <- lme(
  b              ~ section_coverage + age,
  random         = ~ 1 | tank,
  weights        = varIdent(form = ~ 1 | section_coverage),
  data           = df_mod,
  method         = "REML",
  control        = ctrl
)

weighted_chroma <- lme(
  chroma         ~ section_coverage + age,
  random         = ~ 1 | tank,
  weights        = varIdent(form = ~ 1 | section_coverage),
  data           = df_mod,
  method         = "REML",
  control        = ctrl
)

# Sanity check — should see two residual SDs (one per coverage level)

for (lst in list(list(weighted_b, "b*"), list(weighted_L, "L*"),
                 list(weighted_a, "a*"), list(weighted_chroma, "chroma"))) {
  mod  <- lst[[1]]; resp <- lst[[2]]
  sigs <- coef(mod$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE) *
          sigma(mod)
  cat(sprintf("\n%s  — residual SD:  single = %.4f  |  double = %.4f\n",
              resp, sigs["single"], sigs["double"]))
}

## MODEL COMPARISON: AIC / BIC / LOG-LIK

# lmer models refit with ML (not REML) for a fair comparison

# Refit lmer with ML for fair AIC comparison
lmer_L_ml  <- update(full_L,      REML = FALSE)
lmer_a_ml  <- update(full_a,      REML = FALSE)
lmer_b_ml  <- update(full_b,      REML = FALSE)
lmer_chr_ml <- update(full_chroma, REML = FALSE)

# Refit lme with ML
lme_L_ml  <- update(weighted_L,      method = "ML")
lme_a_ml  <- update(weighted_a,      method = "ML")
lme_b_ml  <- update(weighted_b,      method = "ML")
lme_chr_ml <- update(weighted_chroma, method = "ML")

make_comparison_row <- function(resp, homo_ml, weighted_ml) {
  data.frame(
    response     = resp,
    model        = c("Homoscedastic (lmer)", "Weighted varIdent (lme)"),
    df           = c(attr(logLik(homo_ml),     "df"),
                     attr(logLik(weighted_ml), "df")),
    logLik       = round(c(as.numeric(logLik(homo_ml)),
                           as.numeric(logLik(weighted_ml))), 2),
    AIC          = round(c(AIC(homo_ml), AIC(weighted_ml)), 2),
    BIC          = round(c(BIC(homo_ml), BIC(weighted_ml)), 2),
    delta_AIC    = round(c(0, AIC(weighted_ml) - AIC(homo_ml)), 2)
  )
}

comparison_table <- bind_rows(
  make_comparison_row("L*",     lmer_L_ml,   lme_L_ml),
  make_comparison_row("a*",     lmer_a_ml,   lme_a_ml),
  make_comparison_row("b*",     lmer_b_ml,   lme_b_ml),
  make_comparison_row("chroma", lmer_chr_ml, lme_chr_ml)
)

print(comparison_table, row.names = FALSE)

## LEVENE TEST COMPARISON — before vs after weighting

# Create comparison function
levene_compare <- function(resp_label, homo_mod, wtd_mod, data = df_mod) {

  resid_homo <- residuals(homo_mod, type = "response")
  resid_wtd  <- residuals(wtd_mod,  type = "normalized")

  run_levene <- function(resids, grp_var = "section_coverage") {
    df_t <- data |> mutate(resid_abs = abs(resids))
    fit  <- lm(as.formula(paste("resid_abs ~", grp_var)), data = df_t)
    fst  <- summary(fit)$fstatistic
    p    <- pf(fst["value"], fst["numdf"], fst["dendf"], lower.tail = FALSE)
    c(F_val = unname(round(fst["value"], 3)),
      p_val = unname(round(p, 4)))
  }

  homo_lev <- run_levene(resid_homo)
  wtd_lev  <- run_levene(resid_wtd)

  p_vals <- c(homo_lev["p_val"], wtd_lev["p_val"])

  data.frame(
    response      = resp_label,
    model         = c("Homoscedastic", "Weighted varIdent"),
    residual_type = c("raw", "normalised"),
    Levene_F      = c(homo_lev["F_val"], wtd_lev["F_val"]),
    Levene_p      = p_vals
    )
}

levene_table <- bind_rows(
  levene_compare("Lightness", full_L, weighted_L),
  levene_compare("a*", full_a, weighted_a),
  levene_compare("b*", full_b, weighted_b),
  levene_compare("chroma", full_chroma, weighted_chroma)
)

print(levene_table, row.names = FALSE)

## DIAGNOSTIC PLOTS: weighted model residuals

# Build function
diag_weighted <- function(wtd_mod, label, resp, data = df_mod) {

  fitted_vals  <- fitted(wtd_mod)
  resid_norm   <- residuals(wtd_mod, type = "normalized")
  resid_raw    <- residuals(wtd_mod, type = "response")
  sc_cov       <- data$section_coverage

  df_p <- data.frame(
    fitted   = fitted_vals,
    resid    = resid_norm,
    sq_abs   = sqrt(abs(resid_norm)),
    coverage = sc_cov
  )

  p1 <- ggplot(df_p, aes(x = fitted, y = resid, colour = coverage)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    geom_point(alpha = 0.15, size = 0.7) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.6) +
    scale_colour_manual(values = c("single" = "#8DB4C8", "double" = "#2E6E9E"),
                        name = "Coverage") +
    labs(x = "Fitted", y = "Normalised residuals",
         title    = paste0("Residuals vs fitted — ", label, " (weighted)"),
         subtitle = "Colours = coverage groups; spread should now be equal") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8,  colour = "grey40"))

  p2 <- ggplot(df_p, aes(x = fitted, y = sq_abs, colour = coverage)) +
    geom_point(alpha = 0.15, size = 0.7) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 0.8, span = 0.6) +
    scale_colour_manual(values = c("single" = "#8DB4C8", "double" = "#2E6E9E"),
                        name = "Coverage") +
    labs(x = "Fitted", y = expression(sqrt("|Normalised residuals|")),
         title    = paste0("Scale-location — ", label, " (weighted)"),
         subtitle = "Loess per group should now overlap") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8,  colour = "grey40"))

  p3 <- ggplot(df_p, aes(sample = resid)) +
    stat_qq(alpha = 0.2, size = 0.7, colour = "#2E6E9E") +
    stat_qq_line(colour = "#E07B39", linewidth = 0.8) +
    labs(x = "Theoretical quantiles", y = "Normalised residuals",
         title    = paste0("QQ — ", label, " (weighted)"),
         subtitle = "Normality of variance-adjusted residuals") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8,  colour = "grey40"))

  # Box of raw residuals by coverage — direct visual of variance equalisation
  p4 <- ggplot(df_p, aes(x = coverage, y = resid, fill = coverage)) +
    geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4,
                 colour = "grey30", linewidth = 0.4, alpha = 0.8) +
    scale_fill_manual(values = c("single" = "#8DB4C8", "double" = "#2E6E9E"),
                      guide = "none") +
    labs(x = "Section coverage", y = "Normalised residuals",
         title    = paste0("Residuals by coverage — ", label, " (weighted)"),
         subtitle = "Box widths should now be similar") +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8,  colour = "grey40"))

  panel <- (p1 | p2 | p3 | p4) +
    plot_annotation(
      title = paste0("Weighted model diagnostics — ", label),
      theme = theme(plot.title = element_text(size = 13, face = "bold"))
    )

  print(panel)
  ggsave(file.path(fig_dir, sprintf("plot_diagnostics_weighted_%s.png", resp)),
         panel, width = 16, height = 4, dpi = 300)
  cat(sprintf("  Saved: plot_diagnostics_weighted_%s.png\n", resp))
}

# Print diagnostic plots individually
diag_weighted(weighted_L, "L*", "L")
diag_weighted(weighted_a, "a*", "a")
diag_weighted(weighted_b, "b*", "b")
diag_weighted(weighted_chroma, "chroma", "chroma")

# COEFFICIENT COMPARISON: do fixed-effect estimates shift?

# Extract lme co-efficients for plotting
extract_lmer_coefs <- function(mod, response) {
  s <- as.data.frame(coef(summary(mod)))
  data.frame(
    response = response,
    term     = rownames(s),
    estimate = round(s$Estimate,    4),
    se       = round(s$`Std. Error`,4),
    t        = round(s$`t value`,   3),
    p        = round(s$`Pr(>|t|)`,  4),
    model    = "Homoscedastic",
    ci_lo    = round(s$Estimate - 1.96 * s$`Std. Error`, 4),
    ci_hi    = round(s$Estimate + 1.96 * s$`Std. Error`, 4)
  ) |> filter(term != "(Intercept)")
}

extract_lme_coefs <- function(mod, response) {
  s <- as.data.frame(summary(mod)$tTable)
  data.frame(
    response = response,
    term     = rownames(s),
    estimate = round(s$Value,       4),
    se       = round(s$Std.Error,   4),
    t        = round(s$`t-value`,   3),
    p        = round(s$`p-value`,   4),
    model    = "Weighted varIdent",
    ci_lo    = round(s$Value - 1.96 * s$Std.Error, 4),
    ci_hi    = round(s$Value + 1.96 * s$Std.Error, 4)
  ) |> filter(term != "(Intercept)")
}

coef_compare <- bind_rows(
  extract_lmer_coefs(full_L, "Lightness"),
  extract_lme_coefs(weighted_L, "Lightness"),
  extract_lmer_coefs(full_a, "a*"),
  extract_lme_coefs(weighted_a, "a*"),
  extract_lmer_coefs(full_b, "b*"),
  extract_lme_coefs(weighted_b, "b*"),
  extract_lmer_coefs(full_chroma, "chroma"),
  extract_lme_coefs(weighted_chroma, "chroma")
) |>
  mutate(
    response   = factor(response, levels = c("Lightness", "a*", "b*", "chroma")),
    term_clean = case_when(
      term == "section_coveragedouble" ~ "Section: double vs single",
      term == "age31"                   ~ "Age: 31 vs 19 months",
      TRUE                             ~ term
    )
  )

print(coef_compare |>
        dplyr::select(response, term_clean, model, estimate, se, p) |>
        arrange(response, term_clean, model),
      row.names = FALSE)

## Estimates shift slightly but significance does not change - results are robust

# Print final data using weighted model
for (lst in list(list(weighted_L, "Lightness"),
                 list(weighted_a, "a*"),
                 list(weighted_b, "b*"))) {
  mod <- lst[[1]]; response <- lst[[2]]
  
  cat(sprintf("\n %s \n", response))
  
  tt <- summary(mod)$tTable
  print(tt)
  
  ci <- data.frame(
    term  = rownames(tt),
    ci_lo = round(tt[, "Value"] - 1.96 * tt[, "Std.Error"], 4),
    ci_hi = round(tt[, "Value"] + 1.96 * tt[, "Std.Error"], 4)
  )
  cat("95% CI:\n")
  print(ci, row.names = FALSE)
}

### LIKELIHOOD RATIO TESTS MODEL

# Repeat forward selection with weighted model(varIdent)
fit_steps_weighted <- function(response, data = df_mod) {
  ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)
  
  list(
    m0 = lme(as.formula(paste(response, "~ 1")),
             random  = ~ 1 | tank,
             weights = varIdent(form = ~ 1 | section_coverage),
             data    = data, method = "ML", control = ctrl),
    
    m1 = lme(as.formula(paste(response, "~ section_coverage")),
             random  = ~ 1 | tank,
             weights = varIdent(form = ~ 1 | section_coverage),
             data    = data, method = "ML", control = ctrl),
    
    m2 = lme(as.formula(paste(response, "~ age")),
             random  = ~ 1 | tank,
             weights = varIdent(form = ~ 1 | section_coverage),
             data    = data, method = "ML", control = ctrl),
    
    m3 = lme(as.formula(paste(response, "~ section_coverage + age")),
             random  = ~ 1 | tank,
             weights = varIdent(form = ~ 1 | section_coverage),
             data    = data, method = "ML", control = ctrl)
  )
}

# Fit weighted candidate models for each response
models_L_wtd      <- fit_steps_weighted("lightness")
models_a_wtd      <- fit_steps_weighted("a")
models_b_wtd      <- fit_steps_weighted("b")
models_chroma_wtd <- fit_steps_weighted("chroma")

# LRT
lrt_sequential <- function(mods, response) {
  steps <- list(
    list(reduced = mods$m0, full = mods$m1, term_added = "section_coverage"),
    list(reduced = mods$m1, full = mods$m3, term_added = "age")
  )
  lapply(steps, function(s) {
    lt <- anova(s$reduced, s$full)
    
    chi2_val <- lt$`L.Ratio`[2]
    p_val    <- lt$`p-value`[2]
    
    data.frame(
      response   = response,
      term_added = s$term_added,
      model_prev = deparse(formula(s$reduced)),
      chi2       = round(chi2_val, 3),
      df = lt$df[2] - lt$df[1],
      p_value    = round(p_val,    4),
      sig        = case_when(
        p_val < 0.001 ~ "***",
        p_val < 0.01  ~ "**",
        p_val < 0.05  ~ "*",
        p_val < 0.1   ~ ".",
        TRUE          ~ "ns"
      )
    )
  }) |> bind_rows()
}

# Run on weighted model lists
lrt_all_wtd <- bind_rows(
  lrt_sequential(models_L_wtd,      "Lightness"),
  lrt_sequential(models_a_wtd,      "a*"),
  lrt_sequential(models_b_wtd,      "b*"),
  lrt_sequential(models_chroma_wtd, "chroma")
)

for (resp in c("Lightness", "a*", "b*", "chroma")) {
  cat(sprintf("\nLRT (weighted): %s\n", resp))
  tbl <- lrt_all_wtd |>
    dplyr::filter(response == resp) |>
    dplyr::select(term_added, chi2, df, p_value, sig)
  print(tbl, row.names = FALSE)
}

## PRINT WEIGHTED MODEL SUMMARIES

for (lst in list(list(weighted_L, "Lightness"),
                 list(weighted_a, "a*"),
                 list(weighted_b, "b*"),
                 list(weighted_chroma, "chroma"))) {
  mod <- lst[[1]]; response <- lst[[2]]

  cat(sprintf("\n %s \n", response))
  print(summary(mod))

  # Per-group residual SDs
  sigs <- coef(mod$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE) *
          sigma(mod)
  cat(sprintf("  Residual SD:  single = %.4f  |  double = %.4f\n",
              sigs["single"], sigs["double"]))

  r2 <- performance::r2(mod)
  cat(sprintf("  R² marginal:    %.3f\n", r2$R2_marginal))
  cat(sprintf("  R² conditional: %.3f\n", r2$R2_conditional))
  cat("\n")
}

## ESTIMATED MARGINAL MEANS AND PAIRWISE COMPARISONS

for (lst in list(list(weighted_L, "Lightness*"),
                 list(weighted_a, "a*"),
                 list(weighted_b, "b*"),
                 list(weighted_chroma, "chroma"))) {
  mod <- lst[[1]]; response <- lst[[2]]  
  for (term in c("section_coverage", "age")) {
    cat(sprintf("\n %s: %s \n", response, term))
    emm <- emmeans(mod, specs = term)
    print(emm)
    print(pairs(emm, adjust = "tukey"))
  }
}

# Create EMM dataframe for plotting further down

emm_sec_vals <- bind_rows(
  emmeans(weighted_L, specs = "section_coverage") |> as.data.frame() |> mutate(metric = "lightness", level_var = "section_coverage"),
  emmeans(weighted_a, specs = "section_coverage") |> as.data.frame() |> mutate(metric = "a", level_var = "section_coverage"),
  emmeans(weighted_b, specs = "section_coverage") |> as.data.frame() |> mutate(metric = "b", level_var = "section_coverage"),
) |> rename(level = section_coverage)

emm_age_vals <- bind_rows(
  emmeans(weighted_L, specs = "age") |> as.data.frame() |> mutate(metric = "lightness", level_var = "age"),
  emmeans(weighted_a, specs = "age") |> as.data.frame() |> mutate(metric = "a", level_var = "age"),
  emmeans(weighted_b, specs = "age") |> as.data.frame() |> mutate(metric = "b", level_var = "age"),
  emmeans(weighted_chroma, specs = "age") |> as.data.frame() |> mutate(metric = "chroma", level_var = "age")
) |> rename(level = age)

emm_all <- bind_rows(emm_sec_vals, emm_age_vals)

print(emm_all)

# Extract co-efficients for plotting

extract_coefs <- function(mod, response) {
  coefs <- as.data.frame(summary(mod)$tTable)          
  coefs$term <- rownames(coefs)
  coefs |>
    rename(estimate = Value, se = Std.Error,            
           t = `t-value`, p = `p-value`) |>
    filter(term != "(Intercept)") |>
    mutate(
      response   = response,
      ci_lo      = estimate - 1.96 * se,
      ci_hi      = estimate + 1.96 * se,
      sig        = ifelse(p < 0.05, "p < 0.05", "p \u2265 0.05"),
      term_clean = case_when(
        term == "section_coveragedouble" ~ "Section coverage: double vs single",
        term == "age31"                  ~ "Age: 31 vs 19 months",
        TRUE                             ~ term
      ),
      term_group = case_when(
        str_detect(term, "section") ~ "section coverage",
        str_detect(term, "age")     ~ "age",
        TRUE                        ~ "Other"
      )
    )
}

# Build coefficient dataframe from weighted models
coef_df <- bind_rows(
  extract_coefs(weighted_L, "Lightness"),
  extract_coefs(weighted_a, "a* (green - red)"),  
  extract_coefs(weighted_b, "b* (blue - yellow)")     
) |>
  mutate(response = factor(response,
                           levels = c("Lightness",
                                      "a* (green - red)",
                                      "b* (blue - yellow)")))

## FOREST PLOT

sig_cols <- c("p < 0.05" = "#1F618D", "p \u2265 0.05" = "#AED6F1")

p_forest <- ggplot(coef_df,
                   aes(x = estimate,
                       y = reorder(term_clean, estimate),
                       colour = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.25, linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 4) +
  facet_wrap(~response, scales = "free_x", nrow = 1) +
  scale_colour_manual(values = sig_cols, name = NULL) +
  labs(x = "Estimate (95% CI)", y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major   = element_blank(),
    strip.text         = element_text(size = 12, face = "bold"),
    axis.text.y        = element_text(size = 10),
    axis.ticks         = element_line(colour = "black"),
    axis.ticks.length  = unit(0.2, "cm"),
    axis.line          = element_line(colour = "black"),
    legend.position    = "none"
  )

print(p_forest)
ggsave(file.path(fig_dir, "p_forest.png"), plot = p_forest, dpi = 300, width = 12, height = 7, units = "in")

## EMM plots

resp_levels <- c("Lightness", "a* (green - red)", "b* (blue - yellow)")

emm_data <- function(mod, term, response) {
  emmeans(mod, specs = term) |>
    as.data.frame() |>
    mutate(response = factor(response, levels = resp_levels)) |>
    rename(level = !!sym(term))
}

emm_plot <- function(emm_df, x_lab) {
  ggplot(emm_df,
         aes(x = level, y = emmean,
             ymin = lower.CL, ymax = upper.CL)) +
    geom_errorbar(width = 0.18, linewidth = 0.5,
                  position = position_dodge(0.4),
                  colour = "black") +
    geom_point(size = 2.5, position = position_dodge(0.4),
               colour = "black") +
    facet_wrap(~response, scales = "free_y", nrow = 1) +
    labs(x = x_lab, y = "Estimated marginal mean (95% CI)") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_blank(),
      strip.text        = element_text(size = 11, face = "bold"),
      axis.text.x       = element_text(size = 10),
      axis.ticks        = element_line(colour = "black"),
      axis.ticks.length = unit(0.2, "cm"),
      axis.line         = element_line(colour = "black"),
      legend.position   = "none"
    )
}

# Build the EMM dataframes the plots expect (one row per response × level)
emm_sec <- bind_rows(
  emm_data(weighted_L, "section_coverage", "Lightness"),
  emm_data(weighted_a, "section_coverage", "a* (green - red)"),
  emm_data(weighted_b, "section_coverage", "b* (blue - yellow)")
)

emm_age <- bind_rows(
  emm_data(weighted_L, "age", "Lightness"),
  emm_data(weighted_a, "age", "a* (green - red)"),
  emm_data(weighted_b, "age", "b* (blue - yellow)")
)

p_emm_sec <- emm_plot(emm_sec, "Section coverage") + ggtitle("A)") + 
    theme(plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)))
p_emm_age <- emm_plot(emm_age, "Age (months)") + ggtitle("B)") + 
    theme(plot.title = element_text(face = "bold"),
    axis.title.x = element_text(margin = margin(t = 15)),
    axis.title.y = element_text(margin = margin(r = 15)))

print(p_emm_sec)
print(p_emm_age)

ggsave(file.path(fig_dir, "p_emm_sec.png"), plot = p_emm_sec, dpi = 300, width = 10, height = 6, units = "in") 
ggsave(file.path(fig_dir, "p_emm_age.png"), plot = p_emm_age, dpi = 300, width = 10, height = 6, units = "in") 

## script 06: Now that we have full models - create swatches based on EMM values

# Pivot emm_all from full weighted models so that each row = one scenario with L, a, b 

emm_wide <- emm_all |>
  filter(metric %in% c("lightness", "a", "b")) |>
  dplyr::select(level_var, level, metric, emmean) |>
  pivot_wider(names_from = metric, values_from = emmean) |>
  rename(L = lightness)

# Convert CIELAB to RGB colour space
lab_matrix <- as.matrix(emm_wide[, c("L", "a", "b")])

rgb_matrix <- convert_colour(lab_matrix,
                              from       = "lab",
                              to         = "rgb",
                              white_from = "D65",
                              white_to   = "D65") / 255   # farver returns 0–255

# Clamp to [0, 1] (sRGB gamut)
rgb_matrix <- pmax(pmin(rgb_matrix, 1), 0)

emm_wide <- emm_wide |>
  mutate(
    r   = rgb_matrix[, 1],
    g   = rgb_matrix[, 2],
    b_c = rgb_matrix[, 3],   # renamed to avoid clash with b* column
    hex = rgb(r, g, b_c),
    scenario = case_when(
      level_var == "section_coverage" & level == "single" ~ "Single shade",
      level_var == "section_coverage" & level == "double" ~ "Double shade",
      level_var == "age"              & level == "19"     ~ "Age 19 months",
      level_var == "age"              & level == "31"     ~ "Age 31 months"
    ),
    group = if_else(level_var == "section_coverage",
                    "Section Coverage", "Age"),
    lab_label = sprintf("L=%.1f  a*=%.1f  b*=%.1f", L, a, b),
    text_col = if_else(L > 55, "#2C2C2C", "white")
  )

# Plot 
p_swatches <- ggplot(emm_wide,
                     aes(x = scenario, y = 1, fill = hex)) +
  geom_tile(colour = "white", linewidth = 1.5, 
            height = 0.85, width = 0.92) +          
  geom_text(aes(label = scenario, colour = text_col),
            y = 1.18, fontface = "bold", size = 4.5, vjust = 0) +
  geom_text(aes(label = lab_label, colour = text_col, fontface = "bold"),
            y = 0.82, size = 4.5, vjust = 1, family = "mono") +
  geom_text(aes(label = toupper(hex), colour = text_col),
            y = 1.0, size = 4.5, fontface = "bold") +
  scale_fill_identity() +
  scale_colour_identity() +
  facet_wrap(~group, scales = "free_x") +
  scale_y_continuous(limits = c(0.55, 1.45), expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 0.6)) +  
  labs(
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    panel.grid       = element_blank(),
    strip.text       = element_text(size = 16, face = "bold", colour = "#333333"),
    plot.title       = element_text(size = 16, face = "bold"),
    plot.subtitle    = element_text(size = 16,  colour = "grey40"),
    plot.caption     = element_text(size = 16,  colour = "grey55", face = "italic"),
    panel.spacing    = unit(2, "lines")
  )

print(p_swatches)
ggsave(file.path(fig_dir, "plot_emm_colour_swatches.png"), plot = p_swatches, width = 15, height = 10, dpi = 300, bg = "#F8F6F2")

# Script 06: Size metrics — collinearity with age + LMM comparison (age vs length vs area)

ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

### PREPARE DATA

df <- df |>
  mutate(age_num = as.numeric(as.character(age)))

# Filter out missing morphometric data
df <- df |> filter(!(length == 0))

resp_cols <- c("Lightness" = "#5B8FA8", "a*" = "#7A9E5A", "b*" = "#C4A24A")
age_cols   <- c("19" = "#E07B39", "31" = "#5B8FA8")

## HISTOGRAMS: size metrics coloured by age

make_hist <- function(var, xlab, binwidth = NULL, x_limits = NULL, x_breaks = waiver()) {
  if (is.null(binwidth)) {
    rng      <- range(df[[var]], na.rm = TRUE)
    binwidth <- (rng[2] - rng[1]) / 40
  }
  summ <- df |>
    filter(!is.na(.data[[var]])) |>
    group_by(age) |>
    summarise(med = median(.data[[var]], na.rm = TRUE), .groups = "drop")

  ggplot(df |> filter(!is.na(.data[[var]])),
         aes(x = .data[[var]], fill = age, colour = age)) +
    scale_x_continuous(
      limits = x_limits,
      breaks = x_breaks) +
    geom_histogram(binwidth = binwidth, alpha = 0.55,
                   position = "identity", linewidth = 0.2) +
    geom_density(aes(y = after_stat(density) *
                       nrow(df |> filter(!is.na(.data[[var]]))) * binwidth),
                 linewidth = 0.7, alpha = 0) +
    geom_vline(data = summ, aes(xintercept = med, colour = age),
               linetype = "dashed", linewidth = 0.7) +
    scale_fill_manual(values = age_cols, name = "Age (months)") +
    scale_colour_manual(values = age_cols, name = "Age (months)") +
    labs(x = xlab, y = "Count") +
    theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        strip.text = element_text(size = 11),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        plot.title = element_text(size = 13, face = "bold"),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
}

p_len  <- make_hist("length", "Length (mm)", x_limits = c(20, 180), x_breaks = seq(20, 180, 20)) + ggtitle("A)")
p_wid  <- make_hist("width", "Width (mm)", x_limits = c(20, 140), x_breaks = seq(20, 140, 20)) + ggtitle("B)")
p_area <- make_hist("area", "Area (mm²)", binwidth = 300, x_limits = c(0, 15000), x_breaks = seq(0, 15000, 2500)) + ggtitle("C)")

# Patch plot
p_size_hists <- (p_len / p_wid / p_area)

print(p_size_hists)
ggsave(file.path(fig_dir, "plot_size_histograms_by_yc.png"), p_size_hists, width = 9, height = 10, dpi = 300)

### COLLINEARITY: size metrics × each other × age

# Prepare data
size_df <- df |>
  dplyr::select(age_num, length, area, width) |>
  tidyr::drop_na()
 
# Pearson correlation + p-values
cor_r_size <- cor(size_df, method = "pearson")
cor_p_size <- cor_pmat(size_df)
 
# Print results
print(round(cor_r_size, 3))
print(format(cor_p_size, scientific = TRUE))
 
# Plot
size_correlation_p <- ggcorrplot(cor_r_size,
                p.mat      = cor_p_size,
                method     = "square",
                type       = "lower",
                lab        = TRUE,
                lab_size   = 5,
                sig.level  = 0.05,
                insig      = "blank",
                colors     = c("#C0392B", "white", "#1F618D")) +
        theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background = element_rect(fill = "white", colour = NA))
 
print(size_correlation_p)
ggsave(file.path(fig_dir, "size_correlation_p.png"), size_correlation_p, width = 10, height = 10, dpi = 300)

### MODELLING DATASETS
# Model A: section_coverage + age

df_modA <- df |>
  filter(!is.na(section_coverage), !is.na(age),
         !is.na(b), !is.na(a), !is.na(lightness),
         lightness != 0) |>
  mutate(
    section_coverage = factor(section_coverage, ordered = FALSE),
    age              = factor(age,               ordered = FALSE),
    section_coverage = relevel(section_coverage, ref = "single"),
    age              = relevel(age,               ref = "19")
  )

# Model B&C: section_coverage + length_scaled OR area_scaled
str(df)

df_modBC <- df |>
  filter(!is.na(section_coverage),
         !is.na(b), !is.na(a), !is.na(lightness),
         length != 0, lightness != 0) |>
  mutate(
    section_coverage        = factor(section_coverage, ordered = FALSE),
    section_coverage        = relevel(section_coverage, ref = "single"),
    length_scaled           = scale(length)[, 1],
    area_scaled             = scale(area)[, 1]
  )

head(df_modBC,30)

# Model A: section_coverage + age

full_A_L <- lme(lightness ~ section_coverage + age,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modA, method = "REML", control = ctrl)

full_A_a <- lme(a         ~ section_coverage + age,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modA, method = "REML", control = ctrl)

full_A_b <- lme(b         ~ section_coverage + age,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modA, method = "REML", control = ctrl)

full_A_chroma <- lme(chroma    ~ section_coverage + age,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modA, method = "REML", control = ctrl)

# Model B: section_coverage + length_scaled 

full_B_L <- lme(lightness ~ section_coverage + length_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)

full_B_a <- lme(a         ~ section_coverage + length_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)
full_B_b <- lme(b         ~ section_coverage + length_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)

full_B_chroma <- lme(chroma    ~ section_coverage + length_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)

# Model C: section_coverage + area_scaled

full_C_L <- lme(lightness ~ section_coverage + area_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)

full_C_a <- lme(a         ~ section_coverage + area_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)
full_C_b <- lme(b         ~ section_coverage + area_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)

full_C_chroma <- lme(chroma    ~ section_coverage + area_scaled,
                random  = ~ 1 | tank,
                weights = varIdent(form = ~ 1 | section_coverage),
                data    = df_modBC, method = "REML", control = ctrl)

### PRINT SUMMARIES

# Create function to print summary
print_lme_summary <- function(mod, label) {
  cat(sprintf("\n%s\n--- %s ---\n%s\n", strrep("=", 60), label, strrep("=", 60)))
  print(summary(mod))

  # Per-group residual SDs
  sigs <- coef(mod$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE) *
          sigma(mod)
  cat(sprintf("\n  Residual SD:  single = %.4f  |  double = %.4f\n",
              sigs["single"], sigs["double"]))

  # R²
  r2 <- performance::r2(mod)
  cat(sprintf("  R² marginal:    %.3f\n", r2$R2_marginal))
  cat(sprintf("  R² conditional: %.3f\n", r2$R2_conditional))
}

# Print summary for all models
print_lme_summary(full_A_L,      "Model A — L*     | section_coverage + age")
print_lme_summary(full_A_a,      "Model A — a*     | section_coverage + age")
print_lme_summary(full_A_b,      "Model A — b*     | section_coverage + age")
print_lme_summary(full_A_chroma, "Model A — chroma | section_coverage + age")
print_lme_summary(full_B_L,      "Model B — L*     | section_coverage + length_scaled")
print_lme_summary(full_B_a,      "Model B — a*     | section_coverage + length_scaled")
print_lme_summary(full_B_b,      "Model B — b*     | section_coverage + length_scaled")
print_lme_summary(full_B_chroma, "Model B — chroma | section_coverage + length_scaled")
print_lme_summary(full_C_L,      "Model C — L*     | section_coverage + area_scaled")
print_lme_summary(full_C_a,      "Model C — a*     | section_coverage + area_scaled")
print_lme_summary(full_C_b,      "Model C — b*     | section_coverage + area_scaled")
print_lme_summary(full_C_chroma, "Model C — chroma | section_coverage + area_scaled")

### COMPARISON TABLE

refit_ml <- function(mod) update(mod, method = "ML")

# Create function to make summary table
get_row <- function(mod, response, model_label, data) {
  r2   <- performance::r2(mod)
  sigs <- coef(mod$modelStruct$varStruct, unconstrained = FALSE, allCoef = TRUE) *
          sigma(mod)
  mod_ml <- refit_ml(mod)

  data.frame(
    response        = response,
    model           = model_label,
    n               = nrow(data),
    n_tanks         = n_distinct(data$tank),
    AIC             = round(AIC(mod_ml), 1),
    R2_marginal     = round(r2$R2_marginal,    3),
    R2_conditional  = round(r2$R2_conditional, 3),
    resid_SD_single = round(sigs["single"], 4),
    resid_SD_double = round(sigs["double"], 4)
  )
}

comparison_table <- bind_rows(
  get_row(full_A_L,      "L*",     "A: age",                df_modA),
  get_row(full_B_L,      "L*",     "B: length",             df_modBC),
  get_row(full_C_L,      "L*",     "C: area",               df_modBC),
  get_row(full_A_a,      "a*",     "A: age",                df_modA),
  get_row(full_B_a,      "a*",     "B: length",             df_modBC),
  get_row(full_C_a,      "a*",     "C: area",               df_modBC),
  get_row(full_A_b,      "b*",     "A: age",                df_modA),
  get_row(full_B_b,      "b*",     "B: length",             df_modBC),
  get_row(full_C_b,      "b*",     "C: area",               df_modBC),
  get_row(full_A_chroma, "chroma", "A: age",                df_modA),
  get_row(full_B_chroma, "chroma", "B: length",             df_modBC),
  get_row(full_C_chroma, "chroma", "C: area",               df_modBC),
)

print(comparison_table, row.names = FALSE)

# Script 07: Sensitivity test to see if shade cloth and age effect hold up when removing images with calibration below threshold 

# Map each r2 to its respective CIELAB column and define thresholds
r2_cols <- c("L*" = "L_r2", "a*" = "a_r2", "b*" = "b_r2")  
r2_thresholds <- c(0.95, 0.90, 0.80)

# Script 07: Sensitivity test — does removing low-calibration (low R^2) images
#            change the final model's effect sizes, significance or fit?

ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

# --- Explicit filtered datasets (CIELAB rows kept only where calibration R^2 >= threshold) ---
make_filtered <- function(data, r2col, threshold) {
  v <- suppressWarnings(as.numeric(data[[r2col]]))
  data[!is.na(v) & v >= threshold, ]
}

L_filtered <- make_filtered(df_mod, "L_r2", 0.90)   # e.g. lightness, R^2 >= 0.90
a_filtered <- make_filtered(df_mod, "a_r2", 0.90)
b_filtered <- make_filtered(df_mod, "b_r2", 0.90)
# chroma derives from a* and b*, so require BOTH calibrations to pass
chroma_filtered <- {
  va <- suppressWarnings(as.numeric(df_mod$a_r2))
  vb <- suppressWarnings(as.numeric(df_mod$b_r2))
  df_mod[!is.na(va) & !is.na(vb) & va >= 0.90 & vb >= 0.90, ]
}

# --- Refit the FINAL weighted model on any subset ---
fit_final <- function(response, data) {
  lme(as.formula(paste(response, "~ section_coverage + age")),
      random  = ~ 1 | tank,
      weights = varIdent(form = ~ 1 | section_coverage),
      data    = data, method = "REML", control = ctrl)
}

# --- Sweep thresholds and build the comparison table ---
sens_spec <- list(
  "L*"     = list(response = "lightness", r2 = "L_r2"),
  "a*"     = list(response = "a",         r2 = "a_r2"),
  "b*"     = list(response = "b",         r2 = "b_r2"),
  "chroma" = list(response = "chroma",    r2 = c("a_r2", "b_r2"))
)
r2_thresholds <- c(0, 0.80, 0.90, 0.95)   # 0 = unfiltered baseline

sensitivity_one <- function(label, spec, data = df_mod) {
  do.call(rbind, lapply(r2_thresholds, function(thr) {
    if (thr == 0) {
      keep <- rep(TRUE, nrow(data))
    } else {
      keep <- Reduce(`&`, lapply(spec$r2, function(cc) {
        v <- suppressWarnings(as.numeric(data[[cc]]))
        !is.na(v) & v >= thr
      }))
    }
    d  <- droplevels(data[keep, ])
    m  <- fit_final(spec$response, d)
    tt <- summary(m)$tTable
    r2 <- performance::r2(m)
    data.frame(
      response       = label,
      threshold      = if (thr == 0) "none" else sprintf("R2>=%.2f", thr),
      n              = nrow(d),
      n_dropped      = nrow(data) - nrow(d),
      shade_est      = round(tt["section_coveragedouble", "Value"],   3),
      shade_p        = round(tt["section_coveragedouble", "p-value"], 4),
      age_est        = round(tt["age31", "Value"],   3),
      age_p          = round(tt["age31", "p-value"], 4),
      R2_marginal    = round(r2$R2_marginal,    3),
      R2_conditional = round(r2$R2_conditional, 3)
    )
  }))
}

sensitivity_table <- do.call(rbind, Map(sensitivity_one, names(sens_spec), sens_spec))
rownames(sensitivity_table) <- NULL
print(sensitivity_table, row.names = FALSE)

## Effects and significance hold up after removing images with correction values < 0.95, 0.90, and 0.85

# Script 08: Splitting data by age cohort

## NUMERIC SUMMARY — length and area, by age cohort
df |>
  filter(length != 0, !is.na(length)) |>
  pivot_longer(c(length, area), names_to = "metric", values_to = "value") |>
  group_by(metric, age) |>
  summarise(n      = n(),
            mean   = round(mean(value),   1),
            sd     = round(sd(value),     1),
            median = round(median(value), 1),
            min    = round(min(value),    1),
            max    = round(max(value),    1),
            .groups = "drop") |>
  arrange(metric, age) |>
  as.data.frame() |>
  print(row.names = FALSE)

# Sample size and spread sufficient to run individual LMM

## BUILD PER-COHORT MODELLING DATASETS
# Both size predictors (length and area) are scaled within each cohort

# Filter to age = 19 months
df_19 <- df |>
  filter(age == 19,
         !is.na(section_coverage),
         !is.na(b), !is.na(a), !is.na(lightness),
         length != 0, lightness != 0) |>
  mutate(
    section_coverage = relevel(factor(section_coverage, ordered = FALSE), ref = "single"),
    length_scaled    = scale(length)[, 1],
    area_scaled      = scale(area)[, 1]
  )

# Filter to age = 31 months
df_31 <- df |>
  filter(age == 31,
         !is.na(section_coverage),
         !is.na(b), !is.na(a), !is.na(lightness),
         length != 0, lightness != 0) |>
  mutate(
    section_coverage = relevel(factor(section_coverage, ordered = FALSE), ref = "single"),
    length_scaled    = scale(length)[, 1],
    area_scaled      = scale(area)[, 1]
  )

## SANITY CHECK: Are both levels of shade coverage available to permit use of varIdent?
check_coverage <- function(data, label) {
  cat(sprintf("\n--- %s ---\n", label))
  data |>
    distinct(tank, section_coverage) |>
    count(section_coverage) |>
    as.data.frame() |>
    print(row.names = FALSE)
  invisible(nlevels(droplevels(data$section_coverage)))
}

check_coverage(df_19, "Age 19 months")
check_coverage(df_31, "Age 31 months")

## RUN LMM

# set CNTRL
ctrl <- lmeControl(opt = "optim", maxIter = 200, msMaxIter = 200)

# Create functions for fitting models
fit_metric <- function(response, predictor, data) {
  f <- as.formula(paste(response, "~ section_coverage +", predictor))
  lme(f,
      random  = ~ 1 | tank,
      weights = varIdent(form = ~ 1 | section_coverage),
      data    = data, method = "REML", control = ctrl)
}

fit_cohort <- function(data, predictor) {
  list(
    L      = fit_metric("lightness", predictor, data),
    a      = fit_metric("a",         predictor, data),
    b      = fit_metric("b",         predictor, data),
    chroma = fit_metric("chroma",    predictor, data)
  )
}

# Fit models — both size metrics in both cohorts
mods_19_length <- fit_cohort(df_19, "length_scaled")
mods_31_length <- fit_cohort(df_31, "length_scaled")
mods_19_area   <- fit_cohort(df_19, "area_scaled")
mods_31_area   <- fit_cohort(df_31, "area_scaled")

## MODEL SUMMARIES
summarise_cohort <- function(mods, cohort, size_metric) {
  bind_rows(lapply(names(mods), function(m) {
    mod <- mods[[m]]
    tt  <- summary(mod)$tTable
    r2  <- performance::r2(mod)
    data.frame(
      cohort      = cohort,
      size_metric = size_metric,
      response    = m,
      term        = rownames(tt),
      estimate    = round(tt[, "Value"],     3),
      se          = round(tt[, "Std.Error"], 3),
      p           = round(tt[, "p-value"],   4),
      R2_marg     = round(r2$R2_marginal,    3),
      R2_cond     = round(r2$R2_conditional, 3),
      row.names = NULL
    ) |>
      filter(term != "(Intercept)")
  }))
}

# Create summary tables
summary_table_size <- bind_rows(
  summarise_cohort(mods_19_length, "Age 19", "length"),
  summarise_cohort(mods_31_length, "Age 31", "length"),
  summarise_cohort(mods_19_area,   "Age 19", "area"),
  summarise_cohort(mods_31_area,   "Age 31", "area")
) |>
  mutate(
    term = dplyr::recode(term,
                         "section_coveragedouble" = "Coverage: double vs single"),
    term = ifelse(grepl("_scaled$", term), "Size (per SD)", term),
    sig  = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**",
                     p < 0.05  ~ "*",   p < 0.1  ~ ".", TRUE ~ "ns")
  ) |>
  arrange(size_metric, cohort, response, term)

# Print results
print(summary_table_size, row.names = FALSE)

### END OF SCRIPT ###
