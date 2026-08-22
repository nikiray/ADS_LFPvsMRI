# ================================================================
# SECTION 2.2: Participant-level bootstrap test of the
# formal double-dissociation contrast
#
# Tests whether:
#
#   (FW association with aperiodic activity
#      - FW association with burst activity)
#
# differs from:
#
#   (susceptibility association with aperiodic activity
#      - susceptibility association with burst activity)
#
# The bootstrap resamples PARTICIPANTS, retaining both hemispheres
# together within each resampled cluster.
#
# This script is designed to be run after the Section 2.2 primary-analysis
# script, using the hemisphere-level pre-QC data exported by that script.
#
# Required input:
#   results/section_2_2_snc_primary/section_2_2_hemisphere_long_all_preQC.csv
#
#
# Key choices:
# - Composite scores are constructed within the matched complete-case sample.
# - Both MRI markers are entered simultaneously in each domain model.
# - Bootstrap resampling is performed at participant level, retaining both hemispheres.
# - Percentile 95% bootstrap confidence intervals are used for inference.
# - Bootstrap sign proportions are retained as descriptive measures of directional stability.
# ================================================================

library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)

# ------------------------------------------------
# 1. User settings
# ------------------------------------------------

set.seed(20260726)

n_boot <- 5000

primary_dir <- file.path("results", "section_2_2_snc_primary")
input_file <- file.path(
  primary_dir,
  "section_2_2_hemisphere_long_all_preQC.csv"
)

bootstrap_dir <- file.path(
  primary_dir,
  "section_2_2_double_dissociation_bootstrap"
)

dir.create(
  bootstrap_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

led_covariate <- "rs-led_medoff"
snc_dwi_qc_threshold <- 10
exclude_missing_snc_dwi_qc <- TRUE

z_sample <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(x)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

if (!file.exists(input_file)) {
  stop(
    "Required input file not found:\n",
    input_file,
    "\nRun the Section 2.2 primary-analysis script first."
  )
}

long_all <- readr::read_csv(
  input_file,
  show_col_types = FALSE
)

# ------------------------------------------------
# 2. Define component variables
# ------------------------------------------------

aperiodic_components <- c(
  "ap_4_39_offset",
  "ap_4_39_slope"
)

burst_components <- c(
  "lb_occ",
  "lb_bd_median",
  "lb_burstrate"
)

required_cols <- c(
  "sub_id",
  "side",
  "age",
  led_covariate,
  "snc_percvoxremaining",
  "snc_fw",
  "snc_chi",
  "bilateral_cgm_fw",
  "bilateral_cgm_chi",
  aperiodic_components,
  burst_components
)

missing_cols <- setdiff(required_cols, names(long_all))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns in long_all:\n",
    paste(missing_cols, collapse = "\n")
  )
}

# ------------------------------------------------
# 3. Construct the matched complete-case sample
# ------------------------------------------------
#
# Every included hemisphere must have:
#   - both MRI markers
#   - both global MRI covariates
#   - all five LFP component measures
#   - age and residual OFF-medication LED
#   - QC passing at the primary threshold
#


dd_wide <- long_all %>%
  mutate(
    snc_dwi_qc_pass = case_when(
      is.na(snc_percvoxremaining) &
        exclude_missing_snc_dwi_qc ~ FALSE,

      is.na(snc_percvoxremaining) &
        !exclude_missing_snc_dwi_qc ~ TRUE,

      snc_percvoxremaining >= snc_dwi_qc_threshold ~ TRUE,

      TRUE ~ FALSE
    )
  ) %>%
  filter(snc_dwi_qc_pass) %>%
  select(all_of(required_cols)) %>%
  drop_na()

if (
  nrow(dd_wide) < 20 ||
  dplyr::n_distinct(dd_wide$sub_id) < 10
) {
  stop(
    "Too few matched complete observations.\n",
    "Hemisphere rows: ", nrow(dd_wide), "\n",
    "Participants: ", dplyr::n_distinct(dd_wide$sub_id)
  )
}

# ------------------------------------------------
# 4. Standardise component measures and create
#    electrophysiological domain scores
# ------------------------------------------------
#
# The composites summarise the low-frequency aperiodic and low-beta burst
# domains examined in the formal crossover analysis and are constructed before model fitting.
#
# Aperiodic score:
#   mean of standardised 4-39 Hz offset and slope
#
# Burst score:
#   mean of standardised low-beta occupancy, duration and rate
#
# Each composite is then re-standardised within the matched sample.

dd_wide <- dd_wide %>%
  mutate(
    ap_offset_z = z_sample(ap_4_39_offset),
    ap_slope_z  = z_sample(ap_4_39_slope),

    burst_occ_z = z_sample(lb_occ),

    burst_duration_z =
      z_sample(lb_bd_median),

    burst_rate_z =
      z_sample(lb_burstrate),

    aperiodic_score_raw = rowMeans(
      cbind(
        ap_offset_z,
        ap_slope_z
      ),
      na.rm = FALSE
    ),

    burst_score_raw = rowMeans(
      cbind(
        burst_occ_z,
        burst_duration_z,
        burst_rate_z
      ),
      na.rm = FALSE
    ),

    aperiodic_score_z =
      z_sample(aperiodic_score_raw),

    burst_score_z =
      z_sample(burst_score_raw),

    fw_z =
      z_sample(snc_fw),

    chi_z =
      z_sample(snc_chi),

    age_z =
      z_sample(age),

    led_z =
      z_sample(.data[[led_covariate]]),

    global_fw_z =
      z_sample(bilateral_cgm_fw),

    global_chi_z =
      z_sample(bilateral_cgm_chi)
  ) %>%
  drop_na(
    aperiodic_score_z,
    burst_score_z,
    fw_z,
    chi_z,
    age_z,
    led_z,
    global_fw_z,
    global_chi_z
  )

write_csv(
  dd_wide,
  file.path(
    bootstrap_dir,
    "double_dissociation_matched_sample.csv"
  )
)

cat(
  "\nMatched sample used for bootstrap:\n",
  "Hemisphere rows: ", nrow(dd_wide), "\n",
  "Participants: ",
  dplyr::n_distinct(dd_wide$sub_id),
  "\n",
  sep = ""
)

# ------------------------------------------------
# 5. Function to fit the two domain models
# ------------------------------------------------
#
# Both domain models include FW and susceptibility simultaneously so that the
# four marker-by-domain slopes are estimated within the same matched sample.
#
# REML = FALSE is used to improve comparability across bootstrap
# samples and avoid some instability in repeated refitting.

fit_domain_models <- function(dat) {

  fit_aperiodic <- lmer(
    aperiodic_score_z ~
      fw_z +
      chi_z +
      age_z +
      led_z +
      global_fw_z +
      global_chi_z +
      (1 | sub_id),
    data = dat,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )

  fit_burst <- lmer(
    burst_score_z ~
      fw_z +
      chi_z +
      age_z +
      led_z +
      global_fw_z +
      global_chi_z +
      (1 | sub_id),
    data = dat,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )

  b_ap <- fixef(fit_aperiodic)
  b_burst <- fixef(fit_burst)

  fw_domain_difference <-
    b_ap["fw_z"] -
    b_burst["fw_z"]

  chi_domain_difference <-
    b_ap["chi_z"] -
    b_burst["chi_z"]

  double_dissociation <-
    fw_domain_difference -
    chi_domain_difference

  tibble(
    fw_aperiodic =
      unname(b_ap["fw_z"]),

    fw_burst =
      unname(b_burst["fw_z"]),

    chi_aperiodic =
      unname(b_ap["chi_z"]),

    chi_burst =
      unname(b_burst["chi_z"]),

    fw_domain_difference =
      unname(fw_domain_difference),

    chi_domain_difference =
      unname(chi_domain_difference),

    double_dissociation =
      unname(double_dissociation),

    aperiodic_model_singular =
      isSingular(fit_aperiodic, tol = 1e-4),

    burst_model_singular =
      isSingular(fit_burst, tol = 1e-4)
  )
}

# ------------------------------------------------
# 6. Fit the observed-sample models
# ------------------------------------------------

observed_dd <- fit_domain_models(dd_wide)

cat("\nObserved double-dissociation estimates:\n")
print(observed_dd)

write_csv(
  observed_dd,
  file.path(
    bootstrap_dir,
    "observed_double_dissociation_estimates.csv"
  )
)

# ------------------------------------------------
# 7. Participant-level cluster bootstrap
# ------------------------------------------------
#
# Participants are sampled with replacement.
# Both hemispheres belonging to each sampled participant are retained.
#
# Each resampled occurrence receives a new cluster ID so that repeated
# selections of the same original participant are treated as distinct
# bootstrap clusters.

participant_ids <- unique(dd_wide$sub_id)
n_participants <- length(participant_ids)

bootstrap_once <- function(iteration) {

  sampled_ids <- sample(
    participant_ids,
    size = n_participants,
    replace = TRUE
  )

  boot_dat <- map2_dfr(
    sampled_ids,
    seq_along(sampled_ids),
    function(original_id, bootstrap_index) {

      dd_wide %>%
        filter(sub_id == original_id) %>%
        mutate(
          original_sub_id = sub_id,

          sub_id = paste0(
            "bootstrap_",
            bootstrap_index,
            "_",
            original_id
          )
        )
    }
  )

  tryCatch(
    fit_domain_models(boot_dat) %>%
      mutate(
        iteration = iteration,
        status = "OK"
      ),

    error = function(e) {
      tibble(
        fw_aperiodic = NA_real_,
        fw_burst = NA_real_,
        chi_aperiodic = NA_real_,
        chi_burst = NA_real_,
        fw_domain_difference = NA_real_,
        chi_domain_difference = NA_real_,
        double_dissociation = NA_real_,
        aperiodic_model_singular = NA,
        burst_model_singular = NA,
        iteration = iteration,
        status = paste(
          "ERROR:",
          conditionMessage(e)
        )
      )
    }
  )
}

cat(
  "\nRunning ",
  n_boot,
  " participant-level bootstrap samples...\n",
  sep = ""
)

bootstrap_results <- map_dfr(
  seq_len(n_boot),
  bootstrap_once
)

write_csv(
  bootstrap_results,
  file.path(
    bootstrap_dir,
    "participant_bootstrap_all_iterations.csv"
  )
)

successful_bootstrap <- bootstrap_results %>%
  filter(
    status == "OK",
    !is.na(double_dissociation)
  )

failed_bootstrap <- bootstrap_results %>%
  filter(
    status != "OK" |
      is.na(double_dissociation)
  )

cat(
  "\nSuccessful bootstrap samples: ",
  nrow(successful_bootstrap),
  " of ",
  n_boot,
  "\n",
  sep = ""
)

cat(
  "Failed bootstrap samples: ",
  nrow(failed_bootstrap),
  "\n",
  sep = ""
)

if (nrow(successful_bootstrap) < 0.90 * n_boot) {
  warning(
    "Fewer than 90% of bootstrap samples fitted successfully. ",
    "Inspect the failure log before interpreting the interval."
  )
}

write_csv(
  failed_bootstrap,
  file.path(
    bootstrap_dir,
    "participant_bootstrap_failed_iterations.csv"
  )
)

# ------------------------------------------------
# 8. Bootstrap confidence intervals and directional stability
# ------------------------------------------------
#
# Percentile confidence intervals are reported.
#
# Percentile confidence intervals are used for inference.
# The proportions of bootstrap estimates at or below and at or above zero
# are retained as descriptive measures of directional stability.

bootstrap_summary <- successful_bootstrap %>%
  summarise(
    n_requested =
      n_boot,

    n_successful =
      n(),

    observed_double_dissociation =
      observed_dd$double_dissociation[[1]],

    bootstrap_mean =
      mean(
        double_dissociation,
        na.rm = TRUE
      ),

    bootstrap_median =
      median(
        double_dissociation,
        na.rm = TRUE
      ),

    bootstrap_se =
      sd(
        double_dissociation,
        na.rm = TRUE
      ),

    ci_lower_95 =
      unname(
        quantile(
          double_dissociation,
          probs = 0.025,
          na.rm = TRUE
        )
      ),

    ci_upper_95 =
      unname(
        quantile(
          double_dissociation,
          probs = 0.975,
          na.rm = TRUE
        )
      ),

    proportion_at_or_below_zero =
      mean(
        double_dissociation <= 0,
        na.rm = TRUE
      ),

    proportion_at_or_above_zero =
      mean(
        double_dissociation >= 0,
        na.rm = TRUE
      ),

    proportion_same_direction =
      mean(
        sign(double_dissociation) ==
          sign(
            observed_dd$
              double_dissociation[[1]]
          ),
        na.rm = TRUE
      ),

    proportion_aperiodic_models_singular =
      mean(
        aperiodic_model_singular,
        na.rm = TRUE
      ),

    proportion_burst_models_singular =
      mean(
        burst_model_singular,
        na.rm = TRUE
      )
  ) %>%
  mutate(
    ci_excludes_zero =
      ci_lower_95 > 0 |
      ci_upper_95 < 0
  )

cat("\nBootstrap summary:\n")
print(bootstrap_summary)

write_csv(
  bootstrap_summary,
  file.path(
    bootstrap_dir,
    "participant_bootstrap_summary.csv"
  )
)

# ------------------------------------------------
# 9. Bootstrap summaries for all four slopes and
#     both within-marker domain differences
# ------------------------------------------------

parameter_names <- c(
  "fw_aperiodic",
  "fw_burst",
  "chi_aperiodic",
  "chi_burst",
  "fw_domain_difference",
  "chi_domain_difference",
  "double_dissociation"
)

summarise_bootstrap_parameter <- function(parameter_name) {

  values <- successful_bootstrap[[parameter_name]]
  observed_value <- observed_dd[[parameter_name]][1]

  lower_tail <- mean(
    values <= 0,
    na.rm = TRUE
  )

  upper_tail <- mean(
    values >= 0,
    na.rm = TRUE
  )

  ci <- stats::quantile(
    values,
    probs = c(0.025, 0.975),
    na.rm = TRUE,
    names = FALSE
  )

  tibble(
    parameter = parameter_name,
    observed = observed_value,
    bootstrap_mean = mean(values, na.rm = TRUE),
    bootstrap_se = stats::sd(values, na.rm = TRUE),
    ci_lower_95 = ci[1],
    ci_upper_95 = ci[2],
    proportion_at_or_below_zero = lower_tail,
    proportion_at_or_above_zero = upper_tail,
    ci_excludes_zero = (
      ci[1] > 0 ||
      ci[2] < 0
    )
  )
}

parameter_summary <- purrr::map_dfr(
  parameter_names,
  summarise_bootstrap_parameter
)

cat("\nParameter-level bootstrap summary:\n")
print(parameter_summary, n = Inf)

write_csv(
  parameter_summary,
  file.path(
    bootstrap_dir,
    "participant_bootstrap_parameter_summary.csv"
  )
)

# ------------------------------------------------
# 10. Completion message
# ------------------------------------------------

# Save R/package version information for reproducibility
writeLines(
  capture.output(sessionInfo()),
  file.path(bootstrap_dir, "sessionInfo.txt")
)

cat(
  "\nBootstrap double-dissociation analysis complete.\n",
  "Outputs saved in:\n",
  bootstrap_dir,
  "\n\nKey files:\n",
  " - observed_double_dissociation_estimates.csv\n",
  " - participant_bootstrap_summary.csv\n",
  " - participant_bootstrap_parameter_summary.csv\n",
  " - participant_bootstrap_all_iterations.csv\n",
  " - participant_bootstrap_failed_iterations.csv\n",
  sep = ""
)
