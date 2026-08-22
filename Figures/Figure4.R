# ================================================================
# FIGURE 4: MRI-linked STN features and contralateral bradykinesia
# ================================================================
#
# Run after:
#   5.section_2_5_bradykinesia.R
#
# This script reads the exact hemisphere-level dataset and model results saved
# by Script 5. It does not refit the inferential models,
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------

analysis_dir <- file.path("results", "section_2_5_bradykinesia")
hemisphere_file <- file.path(
  analysis_dir, "section_2_5_hemisphere_level_data.csv"
)
results_file <- file.path(
  analysis_dir, "section_2_5_lfp_bradykinesia_all_models.csv"
)
figure_dir <- file.path("results", "figures", "figure_4")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(hemisphere_file, results_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Run the Section 2.5 analysis script first. Missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# Figure 2 palette, plus a distinct PPN colour.
hemisphere_colours <- c(
  "Left STN" = "#2F6FB0",
  "Right STN" = "#C05A16"
)

signature_colours <- c(
  "SNc free water" = "#2F6FB0",
  "SNc susceptibility" = "#C05A16",
  "PPN cAD" = "#2A8F79"
)

# ------------------------------------------------
# 2. Helpers
# ------------------------------------------------

check_cols <- function(dat, cols, object_name) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop(
      object_name, " is missing expected columns:\n",
      paste(missing, collapse = "\n")
    )
  }
}

z_sample <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

format_p <- function(x) {
  if (is.na(x)) return("NA")
  if (x < 0.001) return(formatC(x, format = "e", digits = 2))
  formatC(x, format = "f", digits = ifelse(x < 0.01, 4, 3))
}

# ------------------------------------------------
# 3. Read and validate Script 5 outputs
# ------------------------------------------------

long_brady <- read_csv(hemisphere_file, show_col_types = FALSE)
feature_results <- read_csv(results_file, show_col_types = FALSE)

led_covariate <- "rs-led_medoff"

check_cols(
  long_brady,
  c(
    "sub_id", "side", "age", led_covariate,
    "brady_contralateral_off", "ap_4_39_offset"
  ),
  "Section 2.5 hemisphere-level data"
)

check_cols(
  feature_results,
  c(
    "fdr_family", "predictor", "predictor_var", "display_order",
    "n_hemi", "n_sub", "beta", "se", "p", "q_fdr",
    "ci_lower_95", "ci_upper_95", "fdr_significant"
  ),
  "Section 2.5 model results"
)

if (nrow(feature_results) != 8 || n_distinct(feature_results$predictor_var) != 8) {
  stop("Expected exactly eight unique Section 2.5 feature results.")
}

if (!setequal(unique(long_brady$side), c("left", "right"))) {
  stop("Hemisphere coding failed: both left and right must be present.")
}

# Verify that the stored q values still correspond to the two prespecified
# four-test BH families.
fdr_check <- feature_results %>%
  group_by(fdr_family) %>%
  mutate(q_recalculated = p.adjust(p, method = "BH")) %>%
  ungroup()

if (any(abs(fdr_check$q_recalculated - fdr_check$q_fdr) > 1e-10, na.rm = TRUE)) {
  stop("Stored Section 2.5 q values do not match the prespecified FDR families.")
}

# ------------------------------------------------
# 4. Panel a plotting data
# ------------------------------------------------

# Use the same complete-case sample and standardisation as the corresponding
# Script 5 model. Both variables are residualised for visualisation only;
# inferential statistics come directly from Script 5.
scatter_dat <- long_brady %>%
  select(
    sub_id, side, age, all_of(led_covariate),
    brady_contralateral_off, ap_4_39_offset
  ) %>%
  rename(led_raw = all_of(led_covariate)) %>%
  drop_na() %>%
  mutate(
    brady_z = z_sample(brady_contralateral_off),
    offset_z = z_sample(ap_4_39_offset),
    age_z = z_sample(age),
    led_z = z_sample(led_raw)
  )

brady_adjustment <- lm(brady_z ~ age_z + led_z, data = scatter_dat)
offset_adjustment <- lm(offset_z ~ age_z + led_z, data = scatter_dat)

scatter_dat <- scatter_dat %>%
  mutate(
    adjusted_brady = resid(brady_adjustment),
    adjusted_offset = resid(offset_adjustment),
    hemisphere = factor(
      side,
      levels = c("left", "right"),
      labels = c("Left STN", "Right STN")
    )
  )

offset_result <- feature_results %>%
  filter(predictor_var == "ap_4_39_offset")

if (nrow(offset_result) != 1) {
  stop("Could not uniquely identify the 4-39 Hz aperiodic-offset result.")
}

if (
  nrow(scatter_dat) != offset_result$n_hemi ||
    n_distinct(scatter_dat$sub_id) != offset_result$n_sub
) {
  stop(
    "Panel a sample does not match the saved Script 5 result: ",
    nrow(scatter_dat), " versus ", offset_result$n_hemi, " hemispheres; ",
    n_distinct(scatter_dat$sub_id), " versus ", offset_result$n_sub,
    " participants."
  )
}

# ------------------------------------------------
# 5. Panel b plotting data
# ------------------------------------------------

signature_key <- tibble::tribble(
  ~predictor_var, ~signature, ~display_label,
  "ap_4_39_offset", "SNc free water", "4-39 Hz aperiodic offset",
  "ap_4_39_slope", "SNc free water", "4-39 Hz aperiodic slope",
  "lb_occ", "SNc susceptibility", "Low-beta percentage burst time",
  "lb_bd_median", "SNc susceptibility", "Low-beta median burst duration",
  "lb_burstrate", "SNc susceptibility", "Low-beta burst rate",
  "low_beta_logsm", "SNc susceptibility", "Low-beta oscillatory spectral mass",
  "ap_40_80_offset", "PPN cAD", "40-80 Hz aperiodic offset",
  "ap_40_80_slope", "PPN cAD", "40-80 Hz aperiodic slope"
)

forest_dat <- feature_results %>%
  left_join(signature_key, by = "predictor_var")

if (any(is.na(forest_dat$signature)) || any(is.na(forest_dat$display_label))) {
  stop("At least one Section 2.5 feature could not be assigned to a signature.")
}

forest_dat <- forest_dat %>%
  arrange(display_order) %>%
  mutate(
    display_label = factor(display_label, levels = rev(display_label)),
    significance_label = if_else(
      fdr_significant, "FDR significant", "Not FDR significant"
    )
  )

# ------------------------------------------------
# 6. Figure panels
# ------------------------------------------------

theme_figure <- theme_classic(base_size = 10.5) +
  theme(
    axis.text = element_text(size = 9, colour = "black"),
    axis.title = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
    legend.title = element_text(face = "bold"),
    legend.key.width = grid::unit(1.3, "lines"),
    plot.margin = margin(7, 8, 7, 8)
  )

panel_a <- ggplot(scatter_dat, aes(adjusted_offset, adjusted_brady)) +
  geom_line(
    aes(group = sub_id), colour = "grey55", linewidth = 0.3, alpha = 0.35
  ) +
  geom_point(
    aes(shape = hemisphere, colour = hemisphere),
    size = 2.2, alpha = 0.88
  ) +
  geom_smooth(
    method = "lm", formula = y ~ x, se = TRUE,
    colour = "grey15", fill = "grey72", linewidth = 0.85
  ) +
  scale_shape_manual(values = c("Left STN" = 16, "Right STN" = 17)) +
  scale_colour_manual(values = hemisphere_colours) +
  annotate(
    "label",
    x = Inf, y = -Inf, hjust = 1.05, vjust = -0.45,
    label = paste0(
      "β = ", sprintf("%.3f", offset_result$beta),
      ", p = ", format_p(offset_result$p),
      ", q = ", format_p(offset_result$q_fdr)
    ),
    size = 3.1, label.size = 0.2, fill = "white"
  ) +
  labs(
    title = "Low-frequency aperiodic offset and bradykinesia",
    x = "Adjusted 4-39 Hz aperiodic offset",
    y = "Adjusted contralateral\nOFF-medication bradykinesia",
    shape = "Recording hemisphere",
    colour = "Recording hemisphere"
  ) +
  theme_figure +
  theme(legend.position = "bottom")

panel_b <- ggplot(
  forest_dat,
  aes(
    x = beta, y = display_label,
    xmin = ci_lower_95, xmax = ci_upper_95,
    colour = signature
  )
) +
  geom_vline(
    xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.45
  ) +
  geom_errorbar(
    orientation = "y", height = 0, linewidth = 0.72, alpha = 0.92
  ) +
  geom_point(
    aes(shape = significance_label),
    fill = "white", size = 2.9, stroke = 0.9
  ) +
  scale_colour_manual(values = signature_colours, name = "MRI-linked signature") +
  scale_shape_manual(
    values = c("Not FDR significant" = 21, "FDR significant" = 19),
    name = NULL
  ) +
  labs(
    title = "MRI-linked STN features and bradykinesia",
    x = paste0(
      "Standardised association with contralateral\n",
      "OFF-medication bradykinesia (95% CI)"
    ),
    y = NULL
  ) +
  theme_figure +
  theme(legend.position = "bottom", legend.box = "vertical")

# Give the forest plot slightly more width for its feature labels.
figure_4 <- panel_a + panel_b +
  plot_layout(widths = c(0.92, 1.28), guides = "keep") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13))

# ------------------------------------------------
# 7. Save figure and plotting data
# ------------------------------------------------

ggsave(
  file.path(figure_dir, "Figure_4.png"), figure_4,
  width = 12, height = 5.6, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_4.pdf"), figure_4,
  width = 12, height = 5.6, units = "in", device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_4a.png"), panel_a,
  width = 5.4, height = 4.8, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_4b.png"), panel_b,
  width = 7, height = 5.2, units = "in", dpi = 600, bg = "white"
)

write_csv(scatter_dat, file.path(figure_dir, "Figure_4a_plot_data.csv"))
write_csv(forest_dat, file.path(figure_dir, "Figure_4b_plot_data.csv"))

cat(
  "\nFigure 4 complete.\n",
  "Outputs saved in: ", figure_dir, "\n",
  "Panel a sample: ", nrow(scatter_dat), " hemispheres; ",
  n_distinct(scatter_dat$sub_id), " participants.\n",
  sep = ""
)
