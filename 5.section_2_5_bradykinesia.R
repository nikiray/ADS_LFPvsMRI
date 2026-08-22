# =============================================================================
# SECTION 2.5
# MRI-linked STN electrophysiological features and contralateral
# OFF-medication bradykinesia
#
#
# INPUT FILES
#   data/combined_metrics.csv
#   data/qc_metrics.csv
#
# ANALYSES
#   A. MRI-linked STN LFP features -> contralateral OFF bradykinesia
#   C. Secondary direct MRI -> contralateral OFF bradykinesia
#
# PRIMARY MODEL
#   hemisphere-level linear mixed-effects model:
#
#     z(bradykinesia) ~ z(predictor) + z(age) +
#                       z(residual OFF-medication levodopa exposure) +
#                       (1 | participant)
#
# MRI models additionally include the matched whole-brain MRI covariate where
# appropriate.
#
# STANDARDISATION
#   Outcome, predictor and continuous covariates are z-standardised within the
#   complete-case sample used for each model.
#
# MULTIPLICITY
#   LFP analyses:
#     - four aperiodic features corrected together by BH FDR
#     - three low-beta burst measures plus low-beta oscillatory
#       spectral-mass benchmark corrected together by BH FDR
#
#   Direct MRI analyses:
#     - treated as a separate secondary inferential family
#     - BH FDR applied across the direct MRI tests only
#
# QC
#   The QC threshold is applied when an AD/FW-derived measure is included. 
#
# SINGULAR FITS
#   A corresponding fixed-effects linear model is fitted for every model so
#   that isolated singular random-intercept fits can be checked directly.
#
# =============================================================================


# -----------------------------------------------------------------------------
# 0. Packages
# -----------------------------------------------------------------------------

library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(stringr)
library(lme4)
library(lmerTest)
library(car)


# -----------------------------------------------------------------------------
# 1. User settings
# -----------------------------------------------------------------------------

data_file <- file.path("data", "combined_metrics.csv")
qc_file   <- file.path("data", "qc_metrics.csv")



out_dir <- file.path("results", "section_2_5_bradykinesia")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

led_covariate <- "rs-led_medoff"

expected_cohort_n <- 33L

# SNc voxel-retention QC settings used for AD/FW-derived MRI analyses.
snc_qc_threshold <- 10
exclude_missing_snc_qc <- TRUE


# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

z_sample <- function(x) {
  x <- as.numeric(x)

  if (all(is.na(x))) {
    return(x)
  }

  s <- sd(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }

  as.numeric(
    (x - mean(x, na.rm = TRUE)) / s
  )
}


check_cols <- function(dat, cols) {

  missing <- setdiff(
    cols,
    names(dat)
  )

  if (length(missing) > 0) {
    stop(
      "Missing expected columns:\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }
}


extract_lmer_term <- function(fit, term = "predictor_z") {

  sm <- summary(fit)$coefficients

  if (!term %in% rownames(sm)) {
    return(
      tibble(
        beta = NA_real_,
        se = NA_real_,
        df = NA_real_,
        t = NA_real_,
        p = NA_real_
      )
    )
  }

  tibble(
    beta = unname(sm[term, "Estimate"]),
    se   = unname(sm[term, "Std. Error"]),
    df   = if ("df" %in% colnames(sm)) {
      unname(sm[term, "df"])
    } else {
      NA_real_
    },
    t = unname(sm[term, "t value"]),
    p = unname(sm[term, "Pr(>|t|)"])
  )
}


extract_lm_term <- function(fit, term = "predictor_z") {

  sm <- summary(fit)$coefficients

  if (!term %in% rownames(sm)) {
    return(
      tibble(
        lm_beta = NA_real_,
        lm_se = NA_real_,
        lm_t = NA_real_,
        lm_p = NA_real_
      )
    )
  }

  tibble(
    lm_beta = unname(sm[term, "Estimate"]),
    lm_se   = unname(sm[term, "Std. Error"]),
    lm_t    = unname(sm[term, "t value"]),
    lm_p    = unname(sm[term, "Pr(>|t|)"])
  )
}


get_max_vif <- function(fit) {

  v <- tryCatch(
    car::vif(fit),
    error = function(e) NA
  )

  if (
    length(v) == 1 &&
    is.na(v)
  ) {
    return(NA_real_)
  }

  if (is.matrix(v)) {

    if ("GVIF^(1/(2*Df))" %in% colnames(v)) {
      x <- v[, "GVIF^(1/(2*Df))"]
    } else if ("GVIF" %in% colnames(v)) {
      x <- v[, "GVIF"]
    } else {
      x <- as.numeric(v)
    }

  } else {

    x <- as.numeric(v)
  }

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0) {
    return(NA_real_)
  }

  max(
    x,
    na.rm = TRUE
  )
}


# -----------------------------------------------------------------------------
# 3. Read source data and attach SNc QC
# -----------------------------------------------------------------------------

wide_all <- read_csv(
  data_file,
  show_col_types = FALSE
)


qc_cols <- c(
  "left_snc_percvoxremaining",
  "right_snc_percvoxremaining"
)


if (file.exists(qc_file)) {

  qc <- read_csv(
    qc_file,
    show_col_types = FALSE
  ) %>%
    select(
      sub_id,
      all_of(qc_cols)
    )

  check_cols(
    qc,
    c(
      "sub_id",
      qc_cols
    )
  )

  # qc_metrics.csv is the explicit QC source.
  wide_all <- wide_all %>%
    select(
      -any_of(qc_cols)
    ) %>%
    left_join(
      qc,
      by = "sub_id"
    )

} else {

  stop(
    "Input QC file not found: ",
    qc_file,
    call. = FALSE
  )
}


# -----------------------------------------------------------------------------
# 4. Required raw columns
# -----------------------------------------------------------------------------

lfp_availability_cols <- c(

  "left_stn_4-39-fractal-offset_medsoff",
  "right_stn_4-39-fractal-offset_medsoff",

  "left_stn_4-39-fractal-slope_medsoff",
  "right_stn_4-39-fractal-slope_medsoff",

  "left_stn_40-80-fractal-offset_medsoff",
  "right_stn_40-80-fractal-offset_medsoff",

  "left_stn_40-80-fractal-slope_medsoff",
  "right_stn_40-80-fractal-slope_medsoff",

  "left_stn_13-20-thr-zrobust-occupancy_medoff",
  "right_stn_13-20-thr-zrobust-occupancy_medoff",

  "left_stn_13-20-thr-zrobust-burstrate_medoff",
  "right_stn_13-20-thr-zrobust-burstrate_medoff",

  "left_stn_13-20-thr-zrobust-bd-median_medoff",
  "right_stn_13-20-thr-zrobust-bd-median_medoff",

  "left_stn_13-20-oscillatory-logsm_medoff",
  "right_stn_13-20-oscillatory-logsm_medoff"
)


mri_availability_cols <- c(

  # SNc tissue markers
  "left_snc_fw",
  "right_snc_fw",

  "left_snc_chi",
  "right_snc_chi",

  "left_snc_fwcad",
  "right_snc_fwcad",

  # PPN tissue markers
  "left_ppn_fwcad",
  "right_ppn_fwcad",

  "left_ppn_fw",
  "right_ppn_fw",

  # SNc-striatal structural tract measures
  "left_snc_caud_sw_mu",
  "right_snc_caud_sw_mu",

  "left_snc_puta_sw_mu",
  "right_snc_puta_sw_mu",

  # Matched whole-brain MRI covariates
  "bilateral_cgm_fw",
  "bilateral_cgm_chi",
  "bilateral_cgm_fwcad"
)


required_cols <- c(

  "sub_id",
  "age",
  led_covariate,

  # Lateralised OFF-medication bradykinesia
  "right_updrs-bradykinesia_medoff",
  "left_updrs-bradykinesia_medoff",

  # SNc QC
  qc_cols,

  # LFP and MRI variables
  lfp_availability_cols,
  mri_availability_cols
)


check_cols(
  wide_all,
  required_cols
)


# -----------------------------------------------------------------------------
# 5. Fix participant-level MRI + LFP cohort
# -----------------------------------------------------------------------------
#
# Inclusion is defined before reshaping:
#   participant has at least one MRI value
#   AND at least one STN LFP value.
#
# Individual models then use predictor-specific complete cases within this
# fixed participant cohort.
# -----------------------------------------------------------------------------

cohort_flags <- wide_all %>%
  transmute(

    sub_id,

    n_lfp_nonmissing =
      rowSums(
        !is.na(
          pick(
            all_of(
              lfp_availability_cols
            )
          )
        )
      ),

    n_mri_nonmissing =
      rowSums(
        !is.na(
          pick(
            all_of(
              mri_availability_cols
            )
          )
        )
      ),

    has_any_lfp =
      n_lfp_nonmissing > 0,

    has_any_mri =
      n_mri_nonmissing > 0,

    include_mri_lfp_cohort =
      has_any_lfp &
      has_any_mri
  )


analysis_cohort_ids <- cohort_flags %>%
  filter(
    include_mri_lfp_cohort
  ) %>%
  pull(
    sub_id
  )


if (
  !is.null(expected_cohort_n) &&
  length(analysis_cohort_ids) != expected_cohort_n
) {

  stop(
    paste0(
      "MRI + LFP cohort contains ",
      length(analysis_cohort_ids),
      " participants, not the expected ",
      expected_cohort_n,
      ". Inspect section_2_5_participant_cohort_flags.csv before proceeding."
    ),
    call. = FALSE
  )
}


wide <- wide_all %>%
  filter(
    sub_id %in%
      analysis_cohort_ids
  )


write_csv(
  cohort_flags %>%
    arrange(
      desc(
        include_mri_lfp_cohort
      ),
      sub_id
    ),
  file.path(
    out_dir,
    "section_2_5_participant_cohort_flags.csv"
  )
)


write_csv(
  tibble(
    sub_id =
      sort(
        analysis_cohort_ids
      )
  ),
  file.path(
    out_dir,
    "section_2_5_included_participant_ids.csv"
  )
)


cat(
  "\nFixed MRI + LFP cohort:\n",
  "Participants in input data: ",
  n_distinct(
    wide_all$sub_id
  ),
  "\nParticipants in analysis cohort: ",
  n_distinct(
    wide$sub_id
  ),
  "\n",
  sep = ""
)


# -----------------------------------------------------------------------------
# 6. Construct hemisphere-level dataset
# -----------------------------------------------------------------------------
#
# Clinical lateralisation:
#   left STN / MRI hemisphere  -> right-body bradykinesia
#   right STN / MRI hemisphere -> left-body bradykinesia
# -----------------------------------------------------------------------------

make_side <- function(
  dat,
  side = c(
    "left",
    "right"
  )
) {

  side <- match.arg(
    side
  )


  contralateral_brady_col <-
    if (
      side == "left"
    ) {

      "right_updrs-bradykinesia_medoff"

    } else {

      "left_updrs-bradykinesia_medoff"
    }


  dat %>%
    transmute(

      sub_id =
        sub_id,

      side =
        side,

      # Clinical outcome
      brady_contralateral_off =
        .data[[
contralateral_brady_col
        ]],

      # Participant covariates
      age =
        age,

      led_raw =
        .data[[
led_covariate
        ]],

      # SNc QC
      snc_percvoxremaining =
        .data[[
paste0(
              side,
              "_snc_percvoxremaining"
            )
        ]],

      # ----------------------------------------------------------
      # MRI-linked LFP features
      # ----------------------------------------------------------

      ap_4_39_offset =
        .data[[
paste0(
              side,
              "_stn_4-39-fractal-offset_medsoff"
            )
        ]],

      ap_4_39_slope =
        .data[[
paste0(
              side,
              "_stn_4-39-fractal-slope_medsoff"
            )
        ]],

      lb_occ =
        .data[[
paste0(
              side,
              "_stn_13-20-thr-zrobust-occupancy_medoff"
            )
        ]],

      lb_bd_median =
        .data[[
paste0(
              side,
              "_stn_13-20-thr-zrobust-bd-median_medoff"
            )
        ]],

      lb_burstrate =
        .data[[
paste0(
              side,
              "_stn_13-20-thr-zrobust-burstrate_medoff"
            )
        ]],

      ap_40_80_offset =
        .data[[
paste0(
              side,
              "_stn_40-80-fractal-offset_medsoff"
            )
        ]],

      ap_40_80_slope =
        .data[[
paste0(
              side,
              "_stn_40-80-fractal-slope_medsoff"
            )
        ]],

      # Targeted periodic low-beta benchmark
      low_beta_logsm =
        .data[[
paste0(
              side,
              "_stn_13-20-oscillatory-logsm_medoff"
            )
        ]],

      # ----------------------------------------------------------
      # MRI predictors for secondary direct MRI analyses
      # ----------------------------------------------------------

      snc_fw =
        .data[[
paste0(
              side,
              "_snc_fw"
            )
        ]],

      snc_chi =
        .data[[
paste0(
              side,
              "_snc_chi"
            )
        ]],

      snc_ad =
        .data[[
paste0(
              side,
              "_snc_fwcad"
            )
        ]],

      ppn_ad =
        .data[[
paste0(
              side,
              "_ppn_fwcad"
            )
        ]],

      ppn_fw =
        .data[[
paste0(
              side,
              "_ppn_fw"
            )
        ]],

      snc_caud_sw_mu =
        .data[[
paste0(
              side,
              "_snc_caud_sw_mu"
            )
        ]],

      snc_puta_sw_mu =
        .data[[
paste0(
              side,
              "_snc_puta_sw_mu"
            )
        ]],

      # Whole-brain MRI covariates
      bilateral_cgm_fw =
        bilateral_cgm_fw,

      bilateral_cgm_chi =
        bilateral_cgm_chi,

      bilateral_cgm_fwcad =
        bilateral_cgm_fwcad
    )
}


long_all <- bind_rows(
  make_side(
    wide,
    "left"
  ),
  make_side(
    wide,
    "right"
  )
) %>%
  mutate(

    snc_qc_pass =
      case_when(

        is.na(
          snc_percvoxremaining
        ) &
          exclude_missing_snc_qc ~
          FALSE,

        is.na(
          snc_percvoxremaining
        ) &
          !exclude_missing_snc_qc ~
          TRUE,

        snc_percvoxremaining >=
          snc_qc_threshold ~
          TRUE,

        TRUE ~
          FALSE
      )
  )


write_csv(
  long_all,
  file.path(
    out_dir,
    "section_2_5_hemisphere_level_data.csv"
  )
)


write_csv(
  long_all %>%
    filter(
      !snc_qc_pass
    ) %>%
    select(
      sub_id,
      side,
      snc_percvoxremaining
    ),
  file.path(
    out_dir,
    "section_2_5_snc_qc_excluded_hemispheres.csv"
  )
)


# -----------------------------------------------------------------------------
# 7. Define LFP analyses
# -----------------------------------------------------------------------------
#
# The seven carried-forward features correspond to electrophysiological
# signatures identified by the preceding MRI-LFP analyses.
#
# Low-beta oscillatory spectral mass is retained as the targeted periodic
# benchmark. For the paper-aligned multiplicity analysis it enters the
# four-member low-beta family together with the three burst measures.
# -----------------------------------------------------------------------------

lfp_predictors <- tribble(

  ~analysis_type,
  ~axis,
  ~fdr_family,
  ~predictor_var,
  ~predictor,
  ~display_order,

  "MRI-linked LFP",
  "SNc FW / low-frequency aperiodic",
  "aperiodic",
  "ap_4_39_offset",
  "4-39 Hz aperiodic offset",
  1,

  "MRI-linked LFP",
  "SNc FW / low-frequency aperiodic",
  "aperiodic",
  "ap_4_39_slope",
  "4-39 Hz aperiodic slope",
  2,

  "MRI-linked LFP",
  "SNc susceptibility / low-beta bursts",
  "low_beta",
  "lb_occ",
  "Low-beta percentage burst time",
  3,

  "MRI-linked LFP",
  "SNc susceptibility / low-beta bursts",
  "low_beta",
  "lb_bd_median",
  "Low-beta median burst duration",
  4,

  "MRI-linked LFP",
  "SNc susceptibility / low-beta bursts",
  "low_beta",
  "lb_burstrate",
  "Low-beta burst rate",
  5,

  "MRI-linked LFP",
  "PPN cAD / high-frequency aperiodic",
  "aperiodic",
  "ap_40_80_offset",
  "40-80 Hz aperiodic offset",
  6,

  "MRI-linked LFP",
  "PPN cAD / high-frequency aperiodic",
  "aperiodic",
  "ap_40_80_slope",
  "40-80 Hz aperiodic slope",
  7,

  "Targeted benchmark",
  "Periodic low-beta benchmark",
  "low_beta",
  "low_beta_logsm",
  "Low-beta oscillatory spectral mass",
  8
)


# -----------------------------------------------------------------------------
# 8. Run one LFP -> bradykinesia model
# -----------------------------------------------------------------------------

run_lfp_model <- function(
  analysis_type,
  axis,
  fdr_family,
  predictor_var,
  predictor,
  display_order
) {

  message(
    analysis_type,
    ": ",
    predictor
  )


  model_dat <- long_all %>%
    select(

      sub_id,
      side,

      brady_contralateral_off,

      age,
      led_raw,

      all_of(
        predictor_var
      )
    ) %>%
    rename(

      outcome_raw =
        brady_contralateral_off,

      predictor_raw =
        all_of(
          predictor_var
        )
    ) %>%
    drop_na() %>%
    mutate(

      outcome_z =
        z_sample(
          outcome_raw
        ),

      predictor_z =
        z_sample(
          predictor_raw
        ),

      age_z =
        z_sample(
          age
        ),

      led_z =
        z_sample(
          led_raw
        )
    ) %>%
    drop_na(
      outcome_z,
      predictor_z,
      age_z,
      led_z
    )


  mixed_formula <-
    outcome_z ~
      predictor_z +
      age_z +
      led_z +
      (1 | sub_id)


  fixed_formula <-
    outcome_z ~
      predictor_z +
      age_z +
      led_z


  fit_mixed <- lmer(

    mixed_formula,

    data =
      model_dat,

    REML =
      TRUE,

    control =
      lmerControl(

        optimizer =
          "bobyqa",

        optCtrl =
          list(
            maxfun =
              2e5
          )
      )
  )


  fit_fixed <- lm(

    fixed_formula,

    data =
      model_dat
  )


  mixed <-
    extract_lmer_term(
      fit_mixed
    )


  fixed <-
    extract_lm_term(
      fit_fixed
    )


  tibble(

    analysis_type =
      analysis_type,

    axis =
      axis,

    fdr_family =
      fdr_family,

    predictor =
      predictor,

    predictor_var =
      predictor_var,

    display_order =
      display_order,

    n_hemi =
      nrow(
        model_dat
      ),

    n_sub =
      n_distinct(
        model_dat$sub_id
      ),

    beta =
      mixed$beta,

    se =
      mixed$se,

    df =
      mixed$df,

    t =
      mixed$t,

    p =
      mixed$p,

    singular =
      isSingular(
        fit_mixed,
        tol = 1e-4
      ),

    lm_beta =
      fixed$lm_beta,

    lm_p =
      fixed$lm_p,

    max_vif =
      get_max_vif(
        fit_fixed
      )
  )
}


# -----------------------------------------------------------------------------
# 9. Run LFP models and apply paper-aligned FDR correction
# -----------------------------------------------------------------------------

lfp_results <- pmap_dfr(
  lfp_predictors,
  run_lfp_model
) %>%

  group_by(
    fdr_family
  ) %>%

  mutate(

    q_fdr_family =
      p.adjust(
        p,
        method = "BH"
      )
  ) %>%

  ungroup() %>%

  mutate(

    lower =
      beta -
      1.96 *
      se,

    upper =
      beta +
      1.96 *
      se,

    fdr_significant =
      !is.na(
        q_fdr_family
      ) &
      q_fdr_family <
      0.05
  ) %>%

  arrange(
    display_order
  )


write_csv(
  lfp_results,
  file.path(
    out_dir,
    "section_2_5_lfp_results.csv"
  )
)


write_csv(
  lfp_results %>%
    filter(
      singular
    ),
  file.path(
    out_dir,
    "section_2_5_lfp_singular_model_checks.csv"
  )
)


cat(
  "\nMRI-linked LFP and benchmark results, ",
  "FDR corrected within aperiodic and low-beta families:\n",
  sep = ""
)

print(
  lfp_results,
  n = Inf
)


# -----------------------------------------------------------------------------
# 10. Define secondary direct MRI -> bradykinesia analyses
# -----------------------------------------------------------------------------
#
# These are intentionally secondary and are NOT pooled with the LFP tests.
#
# MRI tissue markers are adjusted for the corresponding bilateral cortical
# grey-matter MRI measure, matching the approach used in the MRI-LFP analyses.
#
# SNc-caudate and SNc-putamen structural tract terms have no matched global
# covariate.
# -----------------------------------------------------------------------------

mri_predictors <- tribble(

  ~axis,
  ~predictor_var,
  ~predictor,
  ~mri_family,
  ~global_var,
  ~apply_snc_qc,

  "SNc FW",
  "snc_fw",
  "SNc free water",
  "SNc_MRI",
  "bilateral_cgm_fw",
  TRUE,

  "SNc susceptibility",
  "snc_chi",
  "SNc susceptibility",
  "SNc_MRI",
  "bilateral_cgm_chi",
  FALSE,

  "SNc cAD",
  "snc_ad",
  "SNc cAD",
  "SNc_MRI",
  "bilateral_cgm_fwcad",
  TRUE,

  "PPN cAD",
  "ppn_ad",
  "PPN cAD",
  "PPN_MRI",
  "bilateral_cgm_fwcad",
  FALSE,

  "PPN FW",
  "ppn_fw",
  "PPN free water",
  "PPN_MRI",
  "bilateral_cgm_fw",
  FALSE,

  "SNc-caudate tract",
  "snc_caud_sw_mu",
  "SNc-caudate structural tract SW mu",
  "SNc_tract",
  NA_character_,
  FALSE,

  "SNc-putamen tract",
  "snc_puta_sw_mu",
  "SNc-putamen structural tract SW mu",
  "SNc_tract",
  NA_character_,
  FALSE
)


# -----------------------------------------------------------------------------
# 11. Run one direct MRI -> bradykinesia model
# -----------------------------------------------------------------------------

run_mri_model <- function(
  axis,
  predictor_var,
  predictor,
  mri_family,
  global_var,
  apply_snc_qc
) {

  message(
    "Secondary direct MRI model: ",
    predictor
  )


  analysis_dat <-
    long_all


  qc_applied <-
    FALSE


  if (
    isTRUE(
      apply_snc_qc
    )
  ) {

    analysis_dat <-
      analysis_dat %>%
      filter(
        snc_qc_pass
      )

    qc_applied <-
      TRUE
  }


  use_global <-
    !is.na(
      global_var
    ) &&
    length(
      global_var
    ) ==
    1 &&
    nzchar(
      global_var
    )


  needed <- c(

    "sub_id",
    "side",

    "brady_contralateral_off",

    "age",
    "led_raw",

    predictor_var
  )


  if (
    use_global
  ) {

    needed <-
      c(
        needed,
        global_var
      )
  }


  check_cols(
    analysis_dat,
    needed
  )


  model_dat <-
    analysis_dat %>%

    select(
      all_of(
        needed
      )
    ) %>%

    mutate(

      outcome_raw =
        brady_contralateral_off,

      predictor_raw =
        .data[[
predictor_var
        ]],

      global_raw =
        if (
          use_global
        ) {

          .data[[
global_var
        ]]

        } else {

          NA_real_
        }
    ) %>%

    drop_na(

      outcome_raw,
      predictor_raw,
      age,
      led_raw,
      sub_id,
      side
    )


  if (
    use_global
  ) {

    model_dat <-
      model_dat %>%
      drop_na(
        global_raw
      )
  }


  model_dat <-
    model_dat %>%
    mutate(

      outcome_z =
        z_sample(
          outcome_raw
        ),

      predictor_z =
        z_sample(
          predictor_raw
        ),

      age_z =
        z_sample(
          age
        ),

      led_z =
        z_sample(
          led_raw
        ),

      global_z =
        if (
          use_global
        ) {

          z_sample(
            global_raw
          )

        } else {

          NA_real_
        }
    )


  if (
    use_global
  ) {

    model_dat <-
      model_dat %>%
      drop_na(

        outcome_z,
        predictor_z,
        age_z,
        led_z,
        global_z
      )

  } else {

    model_dat <-
      model_dat %>%
      drop_na(

        outcome_z,
        predictor_z,
        age_z,
        led_z
      )
  }


  if (
    use_global
  ) {

    mixed_formula <-
      outcome_z ~
        predictor_z +
        age_z +
        led_z +
        global_z +
        (1 | sub_id)

    fixed_formula <-
      outcome_z ~
        predictor_z +
        age_z +
        led_z +
        global_z

  } else {

    mixed_formula <-
      outcome_z ~
        predictor_z +
        age_z +
        led_z +
        (1 | sub_id)

    fixed_formula <-
      outcome_z ~
        predictor_z +
        age_z +
        led_z
  }


  fit_mixed <- lmer(

    mixed_formula,

    data =
      model_dat,

    REML =
      TRUE,

    control =
      lmerControl(

        optimizer =
          "bobyqa",

        optCtrl =
          list(
            maxfun =
              2e5
          )
      )
  )


  fit_fixed <- lm(

    fixed_formula,

    data =
      model_dat
  )


  mixed <-
    extract_lmer_term(
      fit_mixed
    )


  fixed <-
    extract_lm_term(
      fit_fixed
    )


  tibble(

    analysis_type =
      "Secondary direct MRI",

    axis =
      axis,

    mri_family =
      mri_family,

    predictor =
      predictor,

    predictor_var =
      predictor_var,

    whole_brain_covariate =
      if (
        use_global
      ) {

        global_var

      } else {

        "none"
      },

    snc_qc_applied =
      qc_applied,

    n_hemi =
      nrow(
        model_dat
      ),

    n_sub =
      n_distinct(
        model_dat$sub_id
      ),

    beta =
      mixed$beta,

    se =
      mixed$se,

    df =
      mixed$df,

    t =
      mixed$t,

    p =
      mixed$p,

    singular =
      isSingular(
        fit_mixed,
        tol = 1e-4
      ),

    lm_beta =
      fixed$lm_beta,

    lm_p =
      fixed$lm_p,

    max_vif =
      get_max_vif(
        fit_fixed
      ),

    lower =
      mixed$beta -
      1.96 *
      mixed$se,

    upper =
      mixed$beta +
      1.96 *
      mixed$se
  )
}


# -----------------------------------------------------------------------------
# 12. Run secondary MRI analyses and apply separate MRI FDR correction
# -----------------------------------------------------------------------------

mri_results <- pmap_dfr(
  mri_predictors,
  run_mri_model
) %>%

  mutate(

    q_fdr_mri =
      p.adjust(
        p,
        method = "BH"
      ),

    fdr_significant =
      !is.na(
        q_fdr_mri
      ) &
      q_fdr_mri <
      0.05
  ) %>%

  arrange(
    q_fdr_mri,
    p
  )


write_csv(
  mri_results,
  file.path(
    out_dir,
    "section_2_5_secondary_direct_mri_results.csv"
  )
)


write_csv(
  mri_results %>%
    filter(
      singular
    ),
  file.path(
    out_dir,
    "section_2_5_mri_singular_model_checks.csv"
  )
)


cat(
  "\nSecondary direct MRI -> contralateral OFF bradykinesia results:\n"
)

print(
  mri_results,
  n = Inf
)


# -----------------------------------------------------------------------------
# 13. Save combined result table
# -----------------------------------------------------------------------------

combined_results <- bind_rows(

  lfp_results %>%
    transmute(

      analysis_type,
      axis,

      family =
        fdr_family,

      predictor,
      predictor_var,

      whole_brain_covariate =
        "none",

      snc_qc_applied =
        FALSE,

      n_hemi,
      n_sub,
      beta,
      se,
      df,
      t,
      p,

      q =
        q_fdr_family,

      singular,
      lm_beta,
      lm_p,
      max_vif,
      lower,
      upper,
      fdr_significant
    ),

  mri_results %>%
    transmute(

      analysis_type,
      axis,

      family =
        mri_family,

      predictor,
      predictor_var,
      whole_brain_covariate,
      snc_qc_applied,

      n_hemi,
      n_sub,
      beta,
      se,
      df,
      t,
      p,

      q =
        q_fdr_mri,

      singular,
      lm_beta,
      lm_p,
      max_vif,
      lower,
      upper,
      fdr_significant
    )
)


write_csv(
  combined_results,
  file.path(
    out_dir,
    "section_2_5_all_bradykinesia_results.csv"
  )
)


# -----------------------------------------------------------------------------
# 14. Reproducibility record
# -----------------------------------------------------------------------------

writeLines(
  capture.output(
    sessionInfo()
  ),
  file.path(
    out_dir,
    "sessionInfo.txt"
  )
)


cat(
  "\nSection 2.5 bradykinesia analysis complete.\n",
  "Outputs saved in: ",
  out_dir,
  "\n\n",
  "Key output files:\n",
  "  - section_2_5_lfp_results.csv\n",
  "  - section_2_5_secondary_direct_mri_results.csv\n",
  "  - section_2_5_all_bradykinesia_results.csv\n",
  "  - section_2_5_lfp_singular_model_checks.csv\n",
  "  - section_2_5_mri_singular_model_checks.csv\n",
  "  - section_2_5_participant_cohort_flags.csv\n",
  "  - section_2_5_included_participant_ids.csv\n",
  "  - section_2_5_hemisphere_level_data.csv\n",
  "  - section_2_5_snc_qc_excluded_hemispheres.csv\n",
  "  - sessionInfo.txt\n",
  sep = ""
)
