# ================================================================
# FIGURE 3: PPN microstructure and STN aperiodic activity
# ================================================================
#
# Run after:
#   3.section_2_3_ppn_primary.R
#   4.section_2_4_ppn_snc_moderation.R
#
# The script reconstructs plot-level data from the shared input dataset and
# checks the plotted estimates against the saved results from Sections 2.3
# and 2.4. No numerical results are entered by hand.
# ================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(lme4)
  library(lmerTest)
  library(patchwork)
})

# ------------------------------------------------
# 1. Paths and settings
# ------------------------------------------------

data_file <- file.path("data", "combined_metrics.csv")
qc_file <- file.path("data", "qc_metrics.csv")



section_2_3_dir <- file.path("results", "section_2_3_ppn")
section_2_4_dir <- file.path("results", "section_2_4_ppn_snc_moderation")
figure_dir <- file.path("results", "figures", "figure_3")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

primary_file <- file.path(
  section_2_3_dir, "section_2_3_ppn_primary_all_models.csv"
)
interaction_file <- file.path(
  section_2_4_dir, "section_2_4_primary_ppn_interactions.csv"
)

led_covariate <- "rs-led_medoff"
snc_qc_threshold <- 10
exclude_missing_snc_qc <- TRUE
estimate_tolerance <- 1e-6

required_files <- c(data_file, primary_file, interaction_file)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Run the Section 2.3 and 2.4 analysis scripts first. Missing:\n",
    paste(missing_files, collapse = "\n")
  )
}

# ------------------------------------------------
# 2. Helpers
# ------------------------------------------------

z_sample <- function(x) {
  x <- as.numeric(x)
  s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

check_cols <- function(dat, cols) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop("Missing expected columns:\n", paste(missing, collapse = "\n"))
  }
}

extract_term <- function(fit, term = "predictor_z") {
  sm <- summary(fit)$coefficients
  tibble(
    beta = unname(sm[term, "Estimate"]),
    se = unname(sm[term, "Std. Error"]),
    lower = beta - 1.96 * se,
    upper = beta + 1.96 * se
  )
}

fit_mixed <- function(dat, include_snc = FALSE) {
  formula <- if (include_snc) {
    outcome_z ~ predictor_z + control_z + global_z + age_z + led_z +
      snc_z + (1 | sub_id)
  } else {
    outcome_z ~ predictor_z + control_z + global_z + age_z + led_z +
      (1 | sub_id)
  }

  suppressWarnings(
    lmer(
      formula,
      data = dat,
      REML = TRUE,
      control = lmerControl(
        optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)
      )
    )
  )
}

assert_close <- function(observed, expected, label) {
  if (!is.finite(observed) || !is.finite(expected) ||
      abs(observed - expected) > estimate_tolerance) {
    stop(
      label, " does not agree with the saved analysis output.\n",
      "Plotted estimate: ", signif(observed, 8), "\n",
      "Saved estimate: ", signif(expected, 8)
    )
  }
}

# ------------------------------------------------
# 3. Construct the hemisphere-level dataset
# ------------------------------------------------

wide <- read_csv(data_file, show_col_types = FALSE)

if (file.exists(qc_file)) {
  qc <- read_csv(qc_file, show_col_types = FALSE)
  check_cols(
    qc,
    c("sub_id", "left_snc_percvoxremaining", "right_snc_percvoxremaining")
  )
  wide <- wide %>%
    select(-any_of(c(
      "left_snc_percvoxremaining", "right_snc_percvoxremaining"
    ))) %>%
    left_join(
      qc %>% select(
        sub_id, left_snc_percvoxremaining, right_snc_percvoxremaining
      ),
      by = "sub_id"
    )
}

required_wide <- c(
  "sub_id", "age", led_covariate,
  "left_ppn_fwcad", "right_ppn_fwcad",
  "left_ppn_fw", "right_ppn_fw",
  "left_cst_fwcad", "right_cst_fwcad",
  "left_cst_fw", "right_cst_fw",
  "bilateral_cgm_fwcad", "bilateral_cgm_fw",
  "left_snc_fw", "right_snc_fw",
  "left_snc_chi", "right_snc_chi",
  "left_snc_percvoxremaining", "right_snc_percvoxremaining",
  "left_stn_4-39-fractal-slope_medsoff",
  "right_stn_4-39-fractal-slope_medsoff",
  "left_stn_40-80-fractal-slope_medsoff",
  "right_stn_40-80-fractal-slope_medsoff",
  "left_stn_40-80-fractal-offset_medsoff",
  "right_stn_40-80-fractal-offset_medsoff"
)
check_cols(wide, required_wide)

make_side <- function(dat, side = c("left", "right")) {
  side <- match.arg(side)
  p <- paste0(side, "_")
  hemisphere_label <- if (side == "left") "Left" else "Right"

  dat %>% transmute(
    sub_id,
    side,
    # Use a value defined outside transmute(). This avoids both the
    # car::recode()/dplyr::recode() conflict and tidy-evaluation ambiguity.
    hemisphere = hemisphere_label,
    age,
    led = .data[[led_covariate]],
    ppn_ad = .data[[paste0(p, "ppn_fwcad")]],
    ppn_fw = .data[[paste0(p, "ppn_fw")]],
    cst_ad = .data[[paste0(p, "cst_fwcad")]],
    cst_fw = .data[[paste0(p, "cst_fw")]],
    global_ad = bilateral_cgm_fwcad,
    global_fw = bilateral_cgm_fw,
    snc_fw = .data[[paste0(p, "snc_fw")]],
    snc_chi = .data[[paste0(p, "snc_chi")]],
    snc_percvoxremaining =
      .data[[paste0(p, "snc_percvoxremaining")]],
    ap_4_39_slope =
      .data[[paste0(p, "stn_4-39-fractal-slope_medsoff")]],
    ap_40_80_slope =
      .data[[paste0(p, "stn_40-80-fractal-slope_medsoff")]],
    ap_40_80_offset =
      .data[[paste0(p, "stn_40-80-fractal-offset_medsoff")]]
  )
}

long_all <- bind_rows(make_side(wide, "left"), make_side(wide, "right")) %>%
  mutate(
    hemisphere = factor(hemisphere, levels = c("Left", "Right")),
    snc_fw_qc_pass = case_when(
      is.na(snc_percvoxremaining) & exclude_missing_snc_qc ~ FALSE,
      is.na(snc_percvoxremaining) & !exclude_missing_snc_qc ~ TRUE,
      snc_percvoxremaining >= snc_qc_threshold ~ TRUE,
      TRUE ~ FALSE
    )
  )

if (!setequal(as.character(unique(long_all$hemisphere)), c("Left", "Right"))) {
  stop("Hemisphere coding failed: both Left and Right must be present.")
}

# ------------------------------------------------
# 4. Plot data and exact model estimates
# ------------------------------------------------

prepare_model_data <- function(
    outcome_var, predictor_var, control_var, global_var,
    snc_var = NULL, apply_snc_qc = FALSE
) {
  dat <- if (apply_snc_qc) filter(long_all, snc_fw_qc_pass) else long_all
  cols <- c(
    "sub_id", "side", "hemisphere", "age", "led",
    outcome_var, predictor_var, control_var, global_var, snc_var
  )

  dat <- dat %>%
    select(all_of(cols)) %>%
    drop_na() %>%
    rename(
      outcome_raw = all_of(outcome_var),
      predictor_raw = all_of(predictor_var),
      control_raw = all_of(control_var),
      global_raw = all_of(global_var)
    ) %>%
    mutate(
      outcome_z = z_sample(outcome_raw),
      predictor_z = z_sample(predictor_raw),
      control_z = z_sample(control_raw),
      global_z = z_sample(global_raw),
      age_z = z_sample(age),
      led_z = z_sample(led)
    )

  if (!is.null(snc_var)) {
    dat <- dat %>%
      rename(snc_raw = all_of(snc_var)) %>%
      mutate(snc_z = z_sample(snc_raw))
  }
  dat
}

make_adjusted_scatter <- function(outcome_var) {
  dat <- prepare_model_data(
    outcome_var, "ppn_ad", "cst_ad", "global_ad"
  )

  outcome_fit <- lm(
    outcome_z ~ control_z + global_z + age_z + led_z, data = dat
  )
  predictor_fit <- lm(
    predictor_z ~ control_z + global_z + age_z + led_z, data = dat
  )

  dat %>% mutate(
    adjusted_outcome = resid(outcome_fit),
    adjusted_predictor = resid(predictor_fit)
  )
}

coefficient_specs <- tribble(
  ~effect_id, ~effect_label, ~outcome_var, ~predictor_var,
  ~control_var, ~global_var,
  "cad_slope", "PPN cAD: 40-80 Hz slope", "ap_40_80_slope",
  "ppn_ad", "cst_ad", "global_ad",
  "cad_offset", "PPN cAD: 40-80 Hz offset", "ap_40_80_offset",
  "ppn_ad", "cst_ad", "global_ad",
  "fw_slope", "PPN FW: 4-39 Hz slope", "ap_4_39_slope",
  "ppn_fw", "cst_fw", "global_fw"
)

adjustment_specs <- tribble(
  ~model_label, ~snc_var, ~apply_snc_qc,
  "Primary model", NA_character_, FALSE,
  "+ SNc free water", "snc_fw", TRUE,
  "+ SNc susceptibility", "snc_chi", FALSE
)

coefficient_data <- crossing(coefficient_specs, adjustment_specs) %>%
  pmap_dfr(function(
      effect_id, effect_label, outcome_var, predictor_var,
      control_var, global_var, model_label, snc_var, apply_snc_qc
  ) {
    has_snc <- !is.na(snc_var)
    dat <- prepare_model_data(
      outcome_var, predictor_var, control_var, global_var,
      snc_var = if (has_snc) snc_var else NULL,
      apply_snc_qc = apply_snc_qc
    )
    estimate <- extract_term(fit_mixed(dat, include_snc = has_snc))
    bind_cols(
      tibble(
        effect_id, effect_label, outcome_var, predictor_var,
        model_label, n_hemi = nrow(dat), n_sub = n_distinct(dat$sub_id)
      ),
      estimate
    )
  }) %>%
  mutate(
    effect_label = factor(
      effect_label,
      levels = rev(c(
        "PPN cAD: 40-80 Hz slope",
        "PPN cAD: 40-80 Hz offset",
        "PPN FW: 4-39 Hz slope"
      ))
    ),
    model_label = factor(
      model_label,
      levels = c(
        "Primary model", "+ SNc free water", "+ SNc susceptibility"
      )
    )
  )

# Confirm the three primary points match the exact Section 2.3 output.
saved_primary <- read_csv(primary_file, show_col_types = FALSE)
walk(seq_len(nrow(coefficient_specs)), function(i) {
  spec <- coefficient_specs[i, ]
  observed <- coefficient_data %>%
    filter(effect_id == spec$effect_id, model_label == "Primary model") %>%
    pull(beta)
  expected <- saved_primary %>%
    filter(
      predictor_var == spec$predictor_var,
      outcome_var == spec$outcome_var
    ) %>%
    pull(beta)
  if (length(expected) != 1) stop("Could not uniquely match ", spec$effect_id)
  assert_close(observed, expected, paste("Panel c", spec$effect_id))
})

# Moderation model for panel d (same complete-case standardisation as Section 2.4).
interaction_data <- long_all %>%
  select(
    sub_id, side, ppn_ad, snc_chi, cst_ad, global_ad,
    age, led, ap_40_80_slope
  ) %>%
  drop_na() %>%
  mutate(
    outcome_z = z_sample(ap_40_80_slope),
    ppn_z = z_sample(ppn_ad),
    snc_z = z_sample(snc_chi),
    cst_z = z_sample(cst_ad),
    global_z = z_sample(global_ad),
    age_z = z_sample(age),
    led_z = z_sample(led)
  )

interaction_lm <- lm(
  outcome_z ~ ppn_z * snc_z + cst_z + global_z + age_z + led_z,
  data = interaction_data
)

saved_interaction <- read_csv(interaction_file, show_col_types = FALSE) %>%
  filter(
    interaction == "PPN cAD x SNc susceptibility",
    outcome_var == "ap_40_80_slope"
  )
if (nrow(saved_interaction) != 1) {
  stop("Could not uniquely match the panel d interaction result.")
}
assert_close(
  unname(coef(interaction_lm)["ppn_z:snc_z"]),
  saved_interaction$cluster_beta,
  "Panel d interaction"
)

moderator_levels <- c(
  "Lower SNc susceptibility" = -1,
  "Mean SNc susceptibility" = 0,
  "Higher SNc susceptibility" = 1
)

interaction_predictions <- crossing(
  ppn_z = seq(
    min(interaction_data$ppn_z), max(interaction_data$ppn_z), length.out = 150
  ),
  susceptibility_level = names(moderator_levels)
) %>%
  mutate(
    susceptibility_level = factor(
      susceptibility_level, levels = names(moderator_levels)
    ),
    snc_z = unname(moderator_levels[as.character(susceptibility_level)]),
    cst_z = 0, global_z = 0, age_z = 0, led_z = 0
  )

pred <- predict(interaction_lm, newdata = interaction_predictions, se.fit = TRUE)
interaction_predictions <- interaction_predictions %>% mutate(
  fitted = as.numeric(pred$fit),
  se_fit = as.numeric(pred$se.fit),
  lower = fitted - 1.96 * se_fit,
  upper = fitted + 1.96 * se_fit
)

# ------------------------------------------------
# 5. Figure panels
# ------------------------------------------------

theme_figure <- theme_classic(base_size = 10.5) +
  theme(
    plot.title = element_text(face = "bold", size = 11, margin = margin(b = 6)),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    legend.title = element_text(face = "bold"),
    legend.key.width = unit(1.3, "lines"),
    plot.margin = margin(7, 8, 7, 8)
  )

# Match the blue/orange palette used in Figure 2.
hemisphere_colours <- c(
  "Left" = "#2F6FB0",
  "Right" = "#C05A16"
)

scatter_slope <- make_adjusted_scatter("ap_40_80_slope")
scatter_offset <- make_adjusted_scatter("ap_40_80_offset")

make_scatter_panel <- function(dat, title, y_label) {
  ggplot(dat, aes(adjusted_predictor, adjusted_outcome)) +
    geom_line(aes(group = sub_id), colour = "grey55", linewidth = 0.3, alpha = 0.35) +
    geom_point(
      aes(shape = hemisphere, colour = hemisphere),
      size = 2.15, alpha = 0.88
    ) +
    geom_smooth(
      method = "lm", formula = y ~ x, se = TRUE,
      colour = "grey15", fill = "grey72", linewidth = 0.8
    ) +
    scale_shape_manual(values = c(Left = 16, Right = 17)) +
    scale_colour_manual(values = hemisphere_colours) +
    labs(
      title = title, x = "Adjusted PPN cAD", y = y_label,
      shape = "Hemisphere", colour = "Hemisphere"
    ) +
    theme_figure +
    theme(legend.position = "right")
}

panel_a <- make_scatter_panel(
  scatter_slope,
  "PPN cAD and aperiodic slope",
  "Adjusted 40-80 Hz\naperiodic slope"
)

panel_b <- make_scatter_panel(
  scatter_offset,
  "PPN cAD and aperiodic offset",
  "Adjusted 40-80 Hz\naperiodic offset"
)

panel_c <- ggplot(
  coefficient_data,
  aes(
    x = beta, y = effect_label, xmin = lower, xmax = upper,
    shape = model_label
  )
) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.45) +
  geom_errorbar(
    orientation = "y", height = 0,
    position = position_dodge(width = 0.55), linewidth = 0.65
  ) +
  geom_point(position = position_dodge(width = 0.55), size = 2.4) +
  scale_shape_manual(values = c(16, 17, 15)) +
  labs(
    title = "Stability after nigral adjustment",
    x = "Standardised PPN coefficient (95% CI)", y = NULL,
    shape = "Model"
  ) +
  theme_figure +
  theme(legend.position = "bottom")

panel_d <- ggplot() +
  geom_point(
    data = interaction_data,
    aes(ppn_z, outcome_z, colour = snc_z), size = 1.9, alpha = 0.58
  ) +
  geom_ribbon(
    data = interaction_predictions,
    aes(
      x = ppn_z, ymin = lower, ymax = upper,
      fill = susceptibility_level, group = susceptibility_level
    ),
    alpha = 0.12, colour = NA, show.legend = FALSE
  ) +
  geom_line(
    data = interaction_predictions,
    aes(
      x = ppn_z, y = fitted, linetype = susceptibility_level,
      group = susceptibility_level
    ),
    colour = "black", linewidth = 0.9
  ) +
  scale_colour_viridis_c(option = "D", end = 0.9) +
  scale_fill_manual(values = c("#F8766D", "#66C2A5", "#619CFF")) +
  scale_linetype_manual(values = c("dashed", "solid", "dotdash")) +
  labs(
    title = "Moderation by SNc susceptibility",
    x = "PPN cAD (z score)",
    y = "40-80 Hz aperiodic slope (z score)",
    colour = "SNc susceptibility\n(z score)",
    linetype = "Predicted at"
  ) +
  theme_figure +
  theme(legend.position = "right")

figure_3 <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 13))

# ------------------------------------------------
# 6. Save the figure and its plotting data
# ------------------------------------------------

ggsave(
  file.path(figure_dir, "Figure_3.png"), figure_3,
  width = 12, height = 8.6, units = "in", dpi = 600, bg = "white"
)
ggsave(
  file.path(figure_dir, "Figure_3.pdf"), figure_3,
  width = 12, height = 8.6, units = "in", device = cairo_pdf, bg = "white"
)

write_csv(scatter_slope, file.path(figure_dir, "Figure_3a_plot_data.csv"))
write_csv(scatter_offset, file.path(figure_dir, "Figure_3b_plot_data.csv"))
write_csv(coefficient_data, file.path(figure_dir, "Figure_3c_plot_data.csv"))
write_csv(interaction_data, file.path(figure_dir, "Figure_3d_observed_data.csv"))
write_csv(
  interaction_predictions,
  file.path(figure_dir, "Figure_3d_prediction_data.csv")
)

cat(
  "\nFigure 3 complete.\n",
  "Outputs saved in: ", figure_dir, "\n",
  "Panel a/b sample: ", nrow(scatter_slope), " hemispheres; ",
  n_distinct(scatter_slope$sub_id), " participants.\n",
  "Panel d sample: ", nrow(interaction_data), " hemispheres; ",
  n_distinct(interaction_data$sub_id), " participants.\n",
  sep = ""
)
