# ================================================================
# FIGURE 2: Dissociable nigral MRI–STN electrophysiology profiles
#
# Layout:
#   a. Full forest plot of the primary Section 2.2 models
#   b. Adjusted SNc free-water versus 4-39 Hz aperiodic slope
#   c. Adjusted SNc susceptibility versus low-beta burst occupancy
#   d. Formal matched-sample crossover / double-dissociation plot
#
# Assumes the following scripts have already been run in the same R session:
#
#   1. Cleaned Section 2.2 main analysis, creating:
#        results
#        long_all
#        long_fw_qc
#        led_covariate
#        out_dir
#
#   2. Full participant-bootstrap double-dissociation analysis, creating:
#        dd_wide
#        observed_dd
#        parameter_summary
#
# The script saves both PDF and high-resolution PNG versions.
# ================================================================

library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
library(patchwork)
library(scales)

# ------------------------------------------------
# 1. Required objects and packages
# ------------------------------------------------

required_objects <- c(
  "results",
  "dd_wide",
  "observed_dd",
  "parameter_summary",
  "long_fw_qc",
  "long_all",
  "led_covariate",
  "out_dir"
)

missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1))
]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing:\n",
    paste(missing_objects, collapse = "\n"),
    "\n\nRun 1.section_2_2_snc_primary.R and ",
    "2.section_2_2_double_dissociation_bootstrap.R first. The full ",
    "double-dissociation bootstrap script first."
  )
}

figure_dir <- file.path(
  out_dir,
  "figure_2_revised"
)

dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------
# 2. Figure settings
# ------------------------------------------------
#
# Change these colours here if needed.
# They intentionally match the visual identity of the previous figure.

fw_colour <- "#2F6FB0"
chi_colour <- "#C05A16"

fw_pale <- alpha(fw_colour, 0.22)
chi_pale <- alpha(chi_colour, 0.22)

base_size <- 10

theme_figure <- theme_classic(base_size = base_size, base_family = "sans") +
  theme(
    axis.title = element_text(
      colour = "black"
    ),
    axis.text = element_text(
      colour = "black"
    ),
    plot.title = element_text(
      face = "bold",
      size = rel(1.05),
      margin = margin(b = 6)
    ),
    plot.subtitle = element_text(
      size = rel(0.84),
      lineheight = 1.05,
      margin = margin(b = 7)
    ),
    plot.tag = element_text(
      face = "bold",
      size = rel(1.15)
    ),
    legend.title = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    legend.box.just = "left"
  )

# ------------------------------------------------
# 3. Panel a: revised forest plot
# ------------------------------------------------

forest_outcome_order <- c(
  "4-39 Hz aperiodic offset",
  "4-39 Hz aperiodic slope",
  "40-80 Hz aperiodic offset",
  "40-80 Hz aperiodic slope",

  "Low-beta oscillatory log spectral mass",
  "Low-beta percentage burst time",
  "Low-beta burst rate",
  "Low-beta median burst duration",
  "Low-beta median inter-burst interval",

  "High-beta oscillatory log spectral mass",
  "High-beta percentage burst time",
  "High-beta burst rate",
  "High-beta median burst duration",
  "High-beta median inter-burst interval"
)

forest_label_map <- c(
  "4-39 Hz aperiodic offset" =
    "4–39 Hz offset",

  "4-39 Hz aperiodic slope" =
    "4–39 Hz slope",

  "40-80 Hz aperiodic offset" =
    "40–80 Hz offset",

  "40-80 Hz aperiodic slope" =
    "40–80 Hz slope",

  "Low-beta oscillatory log spectral mass" =
    "Low-β spectral mass",

  "Low-beta percentage burst time" =
    "Low-β burst occupancy",

  "Low-beta burst rate" =
    "Low-β burst rate",

  "Low-beta median burst duration" =
    "Low-β burst duration",

  "Low-beta median inter-burst interval" =
    "Low-β inter-burst interval",

  "High-beta oscillatory log spectral mass" =
    "High-β spectral mass",

  "High-beta percentage burst time" =
    "High-β burst occupancy",

  "High-beta burst rate" =
    "High-β burst rate",

  "High-beta median burst duration" =
    "High-β burst duration",

  "High-beta median inter-burst interval" =
    "High-β inter-burst interval"
)

# Updated q_fdr values are used directly, so the final negative
# SNc susceptibility -> low-beta oscillatory spectral-mass association
# is displayed at full opacity when q_fdr < 0.05.

forest_dat <- results %>%
  filter(
    predictor %in% c(
      "SNc free water",
      "SNc susceptibility"
    ),
    outcome %in% forest_outcome_order
  ) %>%
  mutate(
    marker = case_when(
      predictor == "SNc free water" ~
        "SNc free water",

      predictor == "SNc susceptibility" ~
        "SNc susceptibility",

      TRUE ~ predictor
    ),

    outcome_label = unname(
      forest_label_map[outcome]
    ),

    outcome_order = match(
      outcome,
      forest_outcome_order
    ),

    domain = case_when(
      family == "aperiodic" ~
        "Aperiodic activity",

      family == "low_beta" ~
        "Low-β activity",

      family == "high_beta" ~
        "High-β activity",

      TRUE ~ family
    ),

    ci_lower = beta - 1.96 * se,
    ci_upper = beta + 1.96 * se,

    fdr_significant =
      !is.na(q_fdr) &
      q_fdr < 0.05,

    alpha_group = if_else(
      fdr_significant,
      "FDR q < .05",
      "Not FDR-significant"
    ),

    marker = factor(
      marker,
      levels = c(
        "SNc free water",
        "SNc susceptibility"
      )
    )
  ) %>%
  arrange(
    outcome_order,
    marker
  )

if (nrow(forest_dat) == 0) {
  stop(
    "No matching forest-plot results were found. ",
    "Check predictor and outcome labels in `results`."
  )
}

forest_positions <- tibble(
  outcome = forest_outcome_order,
  outcome_label = unname(
    forest_label_map[forest_outcome_order]
  ),
  outcome_order = seq_along(
    forest_outcome_order
  ),
  domain = c(
    rep("Aperiodic activity", 4),
    rep("Low-β activity", 5),
    rep("High-β activity", 5)
  )
)

forest_dat <- forest_dat %>%
  mutate(
    y = max(forest_positions$outcome_order) -
      outcome_order + 1,

    y_dodged = y +
      if_else(
        marker == "SNc free water",
        0.11,
        -0.11
      )
  )

forest_axis <- forest_positions %>%
  mutate(
    y = max(outcome_order) -
      outcome_order + 1
  )

domain_labels <- forest_axis %>%
  group_by(domain) %>%
  summarise(
    y = mean(range(y)),
    .groups = "drop"
  ) %>%
  mutate(
    domain = case_when(
      domain == "Aperiodic activity" ~
        "Aperiodic\nactivity",

      domain == "Low-β activity" ~
        "Low-β\nactivity",

      domain == "High-β activity" ~
        "High-β\nactivity",

      TRUE ~ domain
    )
  )

domain_boundaries <- c(
  mean(c(
    forest_axis$y[
      forest_axis$outcome_order == 4
    ],
    forest_axis$y[
      forest_axis$outcome_order == 5
    ]
  )),
  mean(c(
    forest_axis$y[
      forest_axis$outcome_order == 9
    ],
    forest_axis$y[
      forest_axis$outcome_order == 10
    ]
  ))
)

panel_a_forest <- ggplot(
  forest_dat,
  aes(
    x = beta,
    y = y_dodged,
    colour = marker
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.55,
    colour = "grey35"
  ) +
  geom_hline(
    yintercept = domain_boundaries,
    linewidth = 0.35,
    colour = "grey78"
  ) +
  geom_errorbarh(
    aes(
      xmin = ci_lower,
      xmax = ci_upper,
      alpha = alpha_group
    ),
    height = 0,
    linewidth = 0.75
  ) +
  geom_point(
    aes(
      alpha = alpha_group,
      shape = marker
    ),
    size = 2.6,
    stroke = 0.4
  ) +
  scale_colour_manual(
    values = c(
      "SNc free water" = fw_colour,
      "SNc susceptibility" = chi_colour
    )
  ) +
  scale_shape_manual(
    values = c(
      "SNc free water" = 16,
      "SNc susceptibility" = 17
    )
  ) +
  scale_alpha_manual(
    values = c(
      "FDR q < .05" = 1,
      "Not FDR-significant" = 0.22
    ),
    guide = "none"
  ) +
  scale_y_continuous(
    breaks = forest_axis$y,
    labels = forest_axis$outcome_label,
    expand = expansion(
      mult = c(0.035, 0.035)
    )
  ) +
  labs(
    x = "Standardised association, β (95% CI)",
    y = NULL,
    colour = NULL,
    shape = NULL
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  theme_figure +
  theme(
    legend.position = "top",
    legend.direction = "horizontal",
    axis.text.y = element_text(
      size = rel(0.85)
    ),
    plot.margin = margin(
      t = 4,
      r = 8,
      b = 4,
      l = 8
    )
  )

# Place domain headings in a separate narrow panel so they do not
# compete with the outcome labels for the same margin space.

panel_a_domains <- ggplot(
  domain_labels,
  aes(
    x = 1,
    y = y,
    label = domain
  )
) +
  geom_text(
    fontface = "bold",
    size = 3.15,
    lineheight = 0.9,
    hjust = 0.5
  ) +
  scale_y_continuous(
    limits = range(forest_axis$y) + c(-0.5, 0.5),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    limits = c(0, 2),
    expand = c(0, 0)
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    tag = "a"
  ) +
  theme_void() +
  theme(
    plot.margin = margin(
      t = 4,
      r = 2,
      b = 4,
      l = 2
    )
  )

panel_a_forest <- panel_a_forest +
  theme(
    plot.margin = margin(
      t = 4,
      r = 8,
      b = 4,
      l = 8
    )
  )

panel_a <- patchwork::wrap_plots(
  panel_a_domains,
  panel_a_forest,
  ncol = 2,
  widths = c(0.13, 1.35)
)

# ------------------------------------------------
# 4. Panel d: formal crossover plot
# ------------------------------------------------

required_parameters <- c(
  "fw_aperiodic",
  "fw_burst",
  "chi_aperiodic",
  "chi_burst"
)

missing_parameters <- setdiff(
  required_parameters,
  parameter_summary$parameter
)

if (length(missing_parameters) > 0) {
  stop(
    "The following bootstrap parameters are missing:\n",
    paste(missing_parameters, collapse = "\n")
  )
}

crossover_dat <- parameter_summary %>%
  filter(
    parameter %in% required_parameters
  ) %>%
  mutate(
    marker = case_when(
      str_starts(parameter, "fw_") ~
        "SNc free water",

      str_starts(parameter, "chi_") ~
        "SNc susceptibility",

      TRUE ~ NA_character_
    ),

    domain = case_when(
      str_ends(parameter, "_aperiodic") ~
        "Low-frequency aperiodic",

      str_ends(parameter, "_burst") ~
        "Low-β burst organisation",

      TRUE ~ NA_character_
    ),

    marker = factor(
      marker,
      levels = c(
        "SNc free water",
        "SNc susceptibility"
      )
    ),

    domain = factor(
      domain,
      levels = c(
        "Low-frequency aperiodic",
        "Low-β burst organisation"
      )
    )
  )

dd_row <- parameter_summary %>%
  filter(
    parameter ==
      "double_dissociation"
  )

if (nrow(dd_row) != 1) {
  stop(
    "Could not uniquely identify the double-dissociation ",
    "bootstrap summary."
  )
}

if (!"proportion_same_direction" %in% names(dd_row)) {
  stop(
    "The double-dissociation bootstrap summary does not contain ",
    "`proportion_same_direction`. Re-run the final bootstrap script first."
  )
}

dd_annotation <- paste0(
  "Difference-in-differences β = ",
  number(
    dd_row$observed,
    accuracy = 0.001
  ),
  "\nParticipant-bootstrap 95% CI ",
  number(
    dd_row$ci_lower_95,
    accuracy = 0.001
  ),
  " to ",
  number(
    dd_row$ci_upper_95,
    accuracy = 0.01
  ),
  "; ",
  number(
    100 * dd_row$proportion_same_direction,
    accuracy = 0.1
  ),
  "% retained the observed direction"
)

panel_d <- ggplot(
  crossover_dat,
  aes(
    x = domain,
    y = observed,
    group = marker,
    colour = marker,
    shape = marker
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.45,
    colour = "grey55"
  ) +
  geom_line(
    linewidth = 1.05
  ) +
  geom_errorbar(
    aes(
      ymin = ci_lower_95,
      ymax = ci_upper_95
    ),
    width = 0.075,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.1,
    stroke = 0.5
  ) +
  scale_colour_manual(
    values = c(
      "SNc free water" = fw_colour,
      "SNc susceptibility" = chi_colour
    )
  ) +
  scale_shape_manual(
    values = c(
      "SNc free water" = 16,
      "SNc susceptibility" = 17
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.10, 0.12))
  ) +
  labs(
    x = NULL,
    y = "Standardised MRI–LFP association",
    subtitle = dd_annotation,
    colour = NULL,
    shape = NULL,
    tag = "d"
  ) +
  theme_figure +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(
      size = 8.2,
      lineheight = 1.05,
      margin = margin(b = 6)
    ),
    axis.text.x = element_text(
      size = rel(0.87)
    ),
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

# ------------------------------------------------
# 5. Helper for primary-analysis adjusted scatterplots
# ------------------------------------------------
#
# This avoids formula errors caused by non-syntactic column names such as
# "rs-led_medoff". The requested columns are first renamed to simple internal
# names and are then residualised against the same covariates used in the
# corresponding primary model.

make_primary_adjusted_plot_data <- function(
    dat,
    predictor_var,
    outcome_var,
    global_var
) {

  needed <- c(
    "sub_id",
    "side",
    "age",
    led_covariate,
    predictor_var,
    outcome_var,
    global_var
  )

  missing_cols <- setdiff(
    needed,
    names(dat)
  )

  if (length(missing_cols) > 0) {
    stop(
      "Missing variables for primary adjusted plot:\n",
      paste(missing_cols, collapse = "\n")
    )
  }

  plot_dat <- dat %>%
    select(
      all_of(needed)
    ) %>%
    transmute(
      sub_id = sub_id,
      side = side,
      age_raw = age,
      led_raw = .data[[led_covariate]],
      predictor_raw = .data[[predictor_var]],
      outcome_raw = .data[[outcome_var]],
      global_raw = .data[[global_var]]
    ) %>%
    drop_na()

  predictor_fit <- lm(
    predictor_raw ~ age_raw + led_raw + global_raw,
    data = plot_dat
  )

  outcome_fit <- lm(
    outcome_raw ~ age_raw + led_raw + global_raw,
    data = plot_dat
  )

  plot_dat %>%
    mutate(
      predictor_adjusted = residuals(predictor_fit),
      outcome_adjusted = residuals(outcome_fit)
    )
}

get_primary_result_annotation <- function(
    predictor_name,
    outcome_name
) {

  result_row <- results %>%
    filter(
      predictor == predictor_name,
      outcome == outcome_name
    )

  if (nrow(result_row) != 1) {
    stop(
      "Could not uniquely identify the primary result for:\n",
      predictor_name,
      " -> ",
      outcome_name
    )
  }

  paste0(
    "β = ",
    scales::number(
      result_row$beta,
      accuracy = 0.001
    ),
    " (95% CI ",
    scales::number(
      result_row$beta - 1.96 * result_row$se,
      accuracy = 0.001
    ),
    " to ",
    scales::number(
      result_row$beta + 1.96 * result_row$se,
      accuracy = 0.001
    ),
    ")\nFDR q = ",
    ifelse(
      result_row$q_fdr < 0.001,
      "< 0.001",
      scales::number(
        result_row$q_fdr,
        accuracy = 0.0001
      )
    ),
    "; ",
    result_row$n_hemi,
    " hemispheres / ",
    result_row$n_sub,
    " participants"
  )
}

# ------------------------------------------------
# 6. Panel b: primary SNc free-water association
# ------------------------------------------------

fw_plot_dat <- make_primary_adjusted_plot_data(
  dat = long_fw_qc,
  predictor_var = "snc_fw",
  outcome_var = "ap_4_39_slope",
  global_var = "bilateral_cgm_fw"
)

fw_annotation <- get_primary_result_annotation(
  predictor_name = "SNc free water",
  outcome_name = "4-39 Hz aperiodic slope"
)

panel_b <- ggplot(
  fw_plot_dat,
  aes(
    x = predictor_adjusted,
    y = outcome_adjusted
  )
) +
  geom_line(
    aes(group = sub_id),
    colour = "grey80",
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = fw_colour,
    fill = alpha(fw_colour, 0.18),
    linewidth = 1.05
  ) +
  geom_point(
    colour = fw_colour,
    size = 2.05,
    alpha = 0.82
  ) +
  labs(
    x = "SNc free water\ncovariate-adjusted",
    y = "4–39 Hz aperiodic slope\ncovariate-adjusted",
    subtitle = fw_annotation,
    tag = "b"
  ) +
  theme_figure +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(
      size = 8.2,
      lineheight = 1.05,
      margin = margin(b = 6)
    ),
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

# ------------------------------------------------
# 7. Panel c: primary SNc susceptibility association
# ------------------------------------------------

chi_plot_dat <- make_primary_adjusted_plot_data(
  dat = long_all,
  predictor_var = "snc_chi",
  outcome_var = "lb_occ",
  global_var = "bilateral_cgm_chi"
)

chi_annotation <- get_primary_result_annotation(
  predictor_name = "SNc susceptibility",
  outcome_name = "Low-beta percentage burst time"
)

panel_c <- ggplot(
  chi_plot_dat,
  aes(
    x = predictor_adjusted,
    y = outcome_adjusted
  )
) +
  geom_line(
    aes(group = sub_id),
    colour = "grey80",
    linewidth = 0.45,
    alpha = 0.75
  ) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    colour = chi_colour,
    fill = alpha(chi_colour, 0.18),
    linewidth = 1.05
  ) +
  geom_point(
    colour = chi_colour,
    size = 2.05,
    alpha = 0.82
  ) +
  labs(
    x = "SNc susceptibility\ncovariate-adjusted",
    y = "Low-β burst occupancy\ncovariate-adjusted",
    subtitle = chi_annotation,
    tag = "c"
  ) +
  theme_figure +
  theme(
    legend.position = "none",
    plot.subtitle = element_text(
      size = 8.2,
      lineheight = 1.05,
      margin = margin(b = 6)
    ),
    plot.margin = margin(
      t = 4,
      r = 4,
      b = 4,
      l = 4
    )
  )

# ------------------------------------------------
# 8. Assemble full figure
# ------------------------------------------------
#
# Final panel order:
#   a. Full forest plot
#   b. Primary SNc free water–aperiodic slope scatter
#   c. Primary SNc susceptibility–burst occupancy scatter
#   d. Formal crossover / difference-in-differences plot

top_right <- patchwork::wrap_plots(
  panel_b,
  panel_c,
  ncol = 2,
  widths = c(1, 1)
)

right_column <- patchwork::wrap_plots(
  top_right,
  panel_d,
  ncol = 1,
  heights = c(1, 1.08)
)

figure_body <- patchwork::wrap_plots(
  panel_a,
  right_column,
  ncol = 2,
  widths = c(1.5, 1.05)
)

figure_2 <- figure_body +
  plot_annotation(
    title = paste0(
      "Nigral free water and susceptibility show dissociable associations ",
      "with STN electrophysiology"
    ),

    subtitle =
      paste0(
        "Free water preferentially tracked low-frequency aperiodic activity, ",
        "whereas susceptibility preferentially tracked low-β burst organisation."
      ),

    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 15,
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(
        size = 10.5,
        lineheight = 1.05,
        margin = margin(b = 10)
      ),
      plot.tag = element_text(
        face = "bold",
        size = 13
      ),
      plot.margin = margin(
        t = 12,
        r = 12,
        b = 10,
        l = 12
      )
    )
  )

# ------------------------------------------------
# 9. Save outputs
# ------------------------------------------------

pdf_file <- file.path(
  figure_dir,
  "Figure_2_revised_double_dissociation.pdf"
)

png_file <- file.path(
  figure_dir,
  "Figure_2_revised_double_dissociation.png"
)

ggsave(
  filename = pdf_file,
  plot = figure_2,
  width = 17,
  height = 9.6,
  units = "in",
  device = cairo_pdf,
  bg = "white"
)

ggsave(
  filename = png_file,
  plot = figure_2,
  width = 17,
  height = 9.6,
  units = "in",
  dpi = 400,
  bg = "white"
)

# Save plotting datasets for audit/reproducibility.

write_csv(
  forest_dat,
  file.path(
    figure_dir,
    "Figure_2_panel_a_forest_data.csv"
  )
)

write_csv(
  crossover_dat,
  file.path(
    figure_dir,
    "Figure_2_panel_d_crossover_data.csv"
  )
)

write_csv(
  fw_plot_dat,
  file.path(
    figure_dir,
    "Figure_2_panel_b_primary_FW_slope_data.csv"
  )
)

write_csv(
  chi_plot_dat,
  file.path(
    figure_dir,
    "Figure_2_panel_c_primary_CHI_occupancy_data.csv"
  )
)

cat(
  "\nRevised Figure 2 saved to:\n",
  pdf_file,
  "\n",
  png_file,
  "\n",
  sep = ""
)
