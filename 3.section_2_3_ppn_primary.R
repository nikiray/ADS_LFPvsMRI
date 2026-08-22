# ================================================================
# SECTION 2.3: PPN MRI predictors of OFF-medication STN physiology

# ================================================================
#
# Purpose
# -------
# Reproduces the Section 2.3 analyses:
#   1. Primary PPN cAD and PPN free-water models
#   2. Standalone CST negative-control models
#   3. SNc-adjusted regional-specificity follow-ups
#   4. Leave-one-participant-out (LOPO) sensitivity analyses
#
# Primary models
# --------------
# PPN cAD:
#   z_LFP ~ z_PPN_cAD + z_CST_cAD + z_whole_brain_cAD +
#           z_age + z_residual_OFF_LED + (1 | sub_id)
#
# PPN free water:
#   z_LFP ~ z_PPN_FW + z_CST_FW + z_whole_brain_FW +
#           z_age + z_residual_OFF_LED + (1 | sub_id)
#
# Key choices
# -----------
# - Continuous variables are z-standardised within each complete-case model.
# - FDR correction is applied separately within each MRI predictor and
#   electrophysiological family.
# - The 10% threshold is applied for QC.
# - Corresponding lm fits are retained to assess singular mixed-effects models.
# - LOPO refits remove both hemispheres for one participant while retaining
#   standardisation derived from the original complete-case model.
# ================================================================

library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(car)

# ------------------------------------------------
# 1. User settings
# ------------------------------------------------



data_file <- file.path("data", "combined_metrics.csv")
qc_file   <- file.path("data", "qc_metrics.csv")   # optional



out_dir <- file.path("results", "section_2_3_ppn")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

led_covariate <- "rs-led_medoff"

snc_qc_threshold <- 10
exclude_missing_snc_qc <- TRUE
#NB, no hemisphere in any participant has significant loss of voxels in PPN

# ------------------------------------------------
# 2. Helper functions
# ------------------------------------------------

z_sample <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(x)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

check_cols <- function(dat, cols) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) {
    stop("Missing expected columns:\n", paste(missing, collapse = "\n"))
  }
}

extract_lmer_term <- function(fit, term = "predictor_z") {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) {
    return(tibble(beta = NA_real_, se = NA_real_, df = NA_real_,
                  t = NA_real_, p = NA_real_))
  }
  tibble(
    beta = unname(sm[term, "Estimate"]),
    se   = unname(sm[term, "Std. Error"]),
    df   = if ("df" %in% colnames(sm)) unname(sm[term, "df"]) else NA_real_,
    t    = unname(sm[term, "t value"]),
    p    = unname(sm[term, "Pr(>|t|)"])
  )
}

extract_lm_term <- function(fit, term = "predictor_z") {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) {
    return(tibble(lm_beta = NA_real_, lm_se = NA_real_,
                  lm_t = NA_real_, lm_p = NA_real_))
  }
  tibble(
    lm_beta = unname(sm[term, "Estimate"]),
    lm_se   = unname(sm[term, "Std. Error"]),
    lm_t    = unname(sm[term, "t value"]),
    lm_p    = unname(sm[term, "Pr(>|t|)"])
  )
}

get_max_vif <- function(fit) {
  v <- tryCatch(car::vif(fit), error = function(e) NA)
  if (length(v) == 1 && all(is.na(v))) return(NA_real_)
  if (is.matrix(v) && "GVIF^(1/(2*Df))" %in% colnames(v)) {
    v <- v[, "GVIF^(1/(2*Df))"]
  } else if (is.matrix(v)) {
    v <- as.numeric(v)
  }
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (length(v) == 0) NA_real_ else max(v)
}

# ------------------------------------------------
# 3. Read data and construct hemisphere-level data
# ------------------------------------------------

wide <- read_csv(data_file, show_col_types = FALSE)

# If a separate QC file is supplied, use its whole-SNc voxel-retention columns.
if (file.exists(qc_file)) {
  qc <- read_csv(qc_file, show_col_types = FALSE)
  check_cols(qc, c("sub_id", "left_snc_percvoxremaining",
                   "right_snc_percvoxremaining"))
  wide <- wide %>%
    select(-any_of(c("left_snc_percvoxremaining",
                     "right_snc_percvoxremaining"))) %>%
    left_join(
      qc %>% select(sub_id, left_snc_percvoxremaining,
                    right_snc_percvoxremaining),
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
  "left_snc_percvoxremaining", "right_snc_percvoxremaining"
)
check_cols(wide, required_wide)

make_side <- function(dat, side = c("left", "right")) {
  side <- match.arg(side)
  p <- paste0(side, "_")

  dat %>%
    transmute(
      sub_id = sub_id,
      side = side,
      age = age,
      !!led_covariate := .data[[led_covariate]],

      ppn_ad = .data[[paste0(p, "ppn_fwcad")]],
      ppn_fw = .data[[paste0(p, "ppn_fw")]],
      cst_ad = .data[[paste0(p, "cst_fwcad")]],
      cst_fw = .data[[paste0(p, "cst_fw")]],

      snc_fw  = .data[[paste0(p, "snc_fw")]],
      snc_chi = .data[[paste0(p, "snc_chi")]],
      snc_percvoxremaining = .data[[paste0(p, "snc_percvoxremaining")]],

      bilateral_cgm_ad = bilateral_cgm_fwcad,
      bilateral_cgm_fw = bilateral_cgm_fw,

      ap_4_39_offset = .data[[paste0(p, "stn_4-39-fractal-offset_medsoff")]],
      ap_4_39_slope  = .data[[paste0(p, "stn_4-39-fractal-slope_medsoff")]],
      ap_40_80_offset = .data[[paste0(p, "stn_40-80-fractal-offset_medsoff")]],
      ap_40_80_slope  = .data[[paste0(p, "stn_40-80-fractal-slope_medsoff")]],

      osc_13_20_logsm = .data[[paste0(p, "stn_13-20-oscillatory-logsm_medoff")]],
      osc_21_30_logsm = .data[[paste0(p, "stn_21-30-oscillatory-logsm_medoff")]],

      lb_occ       = .data[[paste0(p, "stn_13-20-thr-zrobust-occupancy_medoff")]],
      lb_burstrate = .data[[paste0(p, "stn_13-20-thr-zrobust-burstrate_medoff")]],
      lb_bd_median = .data[[paste0(p, "stn_13-20-thr-zrobust-bd-median_medoff")]],
      lb_ibi_median = .data[[paste0(p, "stn_13-20-thr-zrobust-ibi-median_medoff")]],

      hb_occ       = .data[[paste0(p, "stn_21-30-thr-zrobust-occupancy_medoff")]],
      hb_burstrate = .data[[paste0(p, "stn_21-30-thr-zrobust-burstrate_medoff")]],
      hb_bd_median = .data[[paste0(p, "stn_21-30-thr-zrobust-bd-median_medoff")]],
      hb_ibi_median = .data[[paste0(p, "stn_21-30-thr-zrobust-ibi-median_medoff")]]
    )
}

long_all <- bind_rows(
  make_side(wide, "left"),
  make_side(wide, "right")
) %>%
  mutate(
    snc_fw_qc_pass = case_when(
      is.na(snc_percvoxremaining) & exclude_missing_snc_qc ~ FALSE,
      is.na(snc_percvoxremaining) & !exclude_missing_snc_qc ~ TRUE,
      snc_percvoxremaining >= snc_qc_threshold ~ TRUE,
      TRUE ~ FALSE
    )
  )

long_snc_fw_qc <- long_all %>% filter(snc_fw_qc_pass)

cat(
  "\nHemisphere rows in full dataset:", nrow(long_all),
  "\nParticipants in full dataset:", n_distinct(long_all$sub_id),
  "\nHemisphere rows after SNc FW QC:", nrow(long_snc_fw_qc),
  "\nParticipants after SNc FW QC:", n_distinct(long_snc_fw_qc$sub_id),
  "\n"
)

# ------------------------------------------------
# 4. Define predictors and electrophysiological families
# ------------------------------------------------

ppn_predictors <- tribble(
  ~predictor_var, ~predictor, ~control_var, ~control_covariate, ~global_var, ~global_covariate,
  "ppn_ad", "PPN cAD", "cst_ad", "ipsilateral CST cAD",
  "bilateral_cgm_ad", "whole-brain GM cAD",
  "ppn_fw", "PPN free water", "cst_fw", "ipsilateral CST free water",
  "bilateral_cgm_fw", "whole-brain GM free water"
)

cst_predictors <- tribble(
  ~predictor_var, ~predictor, ~global_var, ~global_covariate,
  "cst_ad", "CST cAD", "bilateral_cgm_ad", "whole-brain GM cAD",
  "cst_fw", "CST free water", "bilateral_cgm_fw", "whole-brain GM free water"
)

lfp_outcomes <- tribble(
  ~outcome_var, ~outcome, ~family,

  "ap_4_39_offset",  "4-39 Hz aperiodic offset", "aperiodic",
  "ap_4_39_slope",   "4-39 Hz aperiodic slope", "aperiodic",
  "ap_40_80_offset", "40-80 Hz aperiodic offset", "aperiodic",
  "ap_40_80_slope",  "40-80 Hz aperiodic slope", "aperiodic",

  "osc_13_20_logsm", "Low-beta oscillatory log spectral mass", "low_beta",
  "lb_occ",           "Low-beta percentage burst time", "low_beta",
  "lb_burstrate",     "Low-beta burst rate", "low_beta",
  "lb_bd_median",     "Low-beta median burst duration", "low_beta",
  "lb_ibi_median",    "Low-beta median inter-burst interval", "low_beta",

  "osc_21_30_logsm", "High-beta oscillatory log spectral mass", "high_beta",
  "hb_occ",           "High-beta percentage burst time", "high_beta",
  "hb_burstrate",     "High-beta burst rate", "high_beta",
  "hb_bd_median",     "High-beta median burst duration", "high_beta",
  "hb_ibi_median",    "High-beta median inter-burst interval", "high_beta"
)

# ------------------------------------------------
# 5. Model preparation and fitting
# ------------------------------------------------

prepare_model_data <- function(dat, outcome_var, predictor_var, global_var,
                               control_var = NULL, snc_var = NULL) {

  cols <- c("sub_id", "side", "age", led_covariate,
            outcome_var, predictor_var, global_var,
            control_var, snc_var)
  cols <- cols[!is.na(cols) & nzchar(cols)]
  check_cols(dat, cols)

  out <- dat %>%
    select(all_of(cols)) %>%
    rename(
      outcome_raw   = all_of(outcome_var),
      predictor_raw = all_of(predictor_var),
      global_raw    = all_of(global_var),
      led_raw       = all_of(led_covariate)
    )

  if (!is.null(control_var)) out <- out %>% rename(control_raw = all_of(control_var))
  if (!is.null(snc_var))     out <- out %>% rename(snc_raw = all_of(snc_var))

  raw_required <- c("outcome_raw", "predictor_raw", "global_raw",
                    "age", "led_raw", "sub_id")
  if (!is.null(control_var)) raw_required <- c(raw_required, "control_raw")
  if (!is.null(snc_var))     raw_required <- c(raw_required, "snc_raw")

  out <- out %>% drop_na(all_of(raw_required))

  out <- out %>%
    mutate(
      outcome_z   = z_sample(outcome_raw),
      predictor_z = z_sample(predictor_raw),
      global_z    = z_sample(global_raw),
      age_z       = z_sample(age),
      led_z       = z_sample(led_raw)
    )

  if (!is.null(control_var)) out <- out %>% mutate(control_z = z_sample(control_raw))
  if (!is.null(snc_var))     out <- out %>% mutate(snc_z = z_sample(snc_raw))

  out
}



fit_prepared_model <- function(dat) {
  terms <- c(
    "predictor_z",
    if ("control_z" %in% names(dat)) "control_z",
    if ("snc_z" %in% names(dat)) "snc_z",
    "age_z", "led_z", "global_z"
  )

  fixed_formula <- as.formula(
    paste("outcome_z ~", paste(terms, collapse = " + "))
  )
  mixed_formula <- as.formula(
    paste("outcome_z ~", paste(c(terms, "(1 | sub_id)"), collapse = " + "))
  )

  fit_mixed <- tryCatch(
    suppressWarnings(
      lmer(
        mixed_formula,
        data = dat,
        REML = TRUE,
        control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
      )
    ),
    error = function(e) e
  )

  fit_fixed <- tryCatch(lm(fixed_formula, data = dat), error = function(e) e)

  list(mixed = fit_mixed, fixed = fit_fixed)
}

run_model <- function(dat, outcome_var, outcome, family,
                      predictor_var, predictor,
                      global_var, global_covariate,
                      control_var = NULL, control_covariate = NA_character_,
                      snc_var = NULL, snc_covariate = NA_character_,
                      snc_qc_filter_applied = FALSE) {

  model_dat <- prepare_model_data(
    dat, outcome_var, predictor_var, global_var,
    control_var = control_var, snc_var = snc_var
  )

  base <- tibble(
    predictor = predictor,
    predictor_var = predictor_var,
    control_covariate = control_covariate,
    control_var = ifelse(is.null(control_var), NA_character_, control_var),
    global_covariate = global_covariate,
    global_var = global_var,
    snc_covariate = snc_covariate,
    snc_var = ifelse(is.null(snc_var), NA_character_, snc_var),
    family = family,
    outcome = outcome,
    outcome_var = outcome_var,
    snc_qc_filter_applied = snc_qc_filter_applied,
    n_hemi = nrow(model_dat),
    n_sub = n_distinct(model_dat$sub_id)
  )

  fits <- fit_prepared_model(model_dat)

  if (inherits(fits$mixed, "error")) {
    return(base %>% mutate(
      beta = NA_real_, se = NA_real_, df = NA_real_, t = NA_real_, p = NA_real_,
      singular = NA, lm_beta = NA_real_, lm_se = NA_real_,
      lm_t = NA_real_, lm_p = NA_real_, max_vif = NA_real_,
      status = paste0("lmer error: ", fits$mixed$message)
    ))
  }

  ms <- extract_lmer_term(fits$mixed)

  if (inherits(fits$fixed, "error")) {
    ls <- tibble(lm_beta = NA_real_, lm_se = NA_real_, lm_t = NA_real_, lm_p = NA_real_)
    vif <- NA_real_
  } else {
    ls <- extract_lm_term(fits$fixed)
    vif <- get_max_vif(fits$fixed)
  }

  base %>% mutate(
    beta = ms$beta, se = ms$se, df = ms$df, t = ms$t, p = ms$p,
    singular = isSingular(fits$mixed, tol = 1e-4),
    lm_beta = ls$lm_beta, lm_se = ls$lm_se, lm_t = ls$lm_t, lm_p = ls$lm_p,
    max_vif = vif,
    status = "OK"
  )
}

# ------------------------------------------------
# 6. Primary PPN models
# ------------------------------------------------

ppn_primary_results <- crossing(ppn_predictors, lfp_outcomes) %>%
  pmap_dfr(function(predictor_var, predictor, control_var, control_covariate,
                    global_var, global_covariate, outcome_var, outcome, family) {
    run_model(
      dat = long_all,
      outcome_var = outcome_var, outcome = outcome, family = family,
      predictor_var = predictor_var, predictor = predictor,
      global_var = global_var, global_covariate = global_covariate,
      control_var = control_var, control_covariate = control_covariate
    )
  }) %>%
  group_by(predictor, family) %>%
  mutate(q_fdr = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(predictor, family, q_fdr, p)

ppn_fdr <- ppn_primary_results %>%
  filter(!is.na(q_fdr), q_fdr < 0.05) %>%
  arrange(q_fdr, p)

write_csv(ppn_primary_results, file.path(out_dir, "section_2_3_ppn_primary_all_models.csv"))
write_csv(ppn_fdr, file.path(out_dir, "section_2_3_ppn_primary_fdr_q_lt_0_05.csv"))

cat("\nPrimary PPN findings surviving FDR correction:\n")
print(
  ppn_fdr %>%
    select(predictor, outcome, n_hemi, n_sub, beta, se, t, p, q_fdr,
           singular, lm_beta, lm_p),
  n = Inf
)

# ------------------------------------------------
# 7. Standalone CST negative-control models
# ------------------------------------------------

cst_results <- crossing(cst_predictors, lfp_outcomes) %>%
  pmap_dfr(function(predictor_var, predictor, global_var, global_covariate,
                    outcome_var, outcome, family) {
    run_model(
      dat = long_all,
      outcome_var = outcome_var, outcome = outcome, family = family,
      predictor_var = predictor_var, predictor = predictor,
      global_var = global_var, global_covariate = global_covariate
    )
  }) %>%
  group_by(predictor, family) %>%
  mutate(q_fdr = p.adjust(p, method = "BH")) %>%
  ungroup() %>%
  arrange(predictor, family, q_fdr, p)

write_csv(cst_results, file.path(out_dir, "section_2_3_cst_negative_control_all_models.csv"))
write_csv(
  cst_results %>% filter(!is.na(q_fdr), q_fdr < 0.05),
  file.path(out_dir, "section_2_3_cst_negative_control_fdr_q_lt_0_05.csv")
)

# ------------------------------------------------
# 8. SNc regional-specificity follow-ups
# ------------------------------------------------

primary_effects <- ppn_fdr %>%
  select(
    predictor, predictor_var, control_covariate, control_var,
    global_covariate, global_var, family, outcome, outcome_var,
    primary_beta = beta, primary_se = se, primary_t = t, primary_p = p,
    primary_q = q_fdr, primary_n_hemi = n_hemi, primary_n_sub = n_sub
  )

snc_specs <- tribble(
  ~snc_var, ~snc_covariate, ~requires_snc_fw_qc,
  "snc_fw",  "ipsilateral SNc free water", TRUE,
  "snc_chi", "ipsilateral SNc susceptibility", FALSE
)

if (nrow(primary_effects) > 0) {

  specificity_results <- crossing(primary_effects, snc_specs) %>%
    pmap_dfr(function(
      predictor, predictor_var, control_covariate, control_var,
      global_covariate, global_var, family, outcome, outcome_var,
      primary_beta, primary_se, primary_t, primary_p, primary_q,
      primary_n_hemi, primary_n_sub,
      snc_var, snc_covariate, requires_snc_fw_qc
    ) {

      analysis_dat <- if (requires_snc_fw_qc) long_snc_fw_qc else long_all

      adjusted_dat <- prepare_model_data(
        analysis_dat, outcome_var, predictor_var, global_var,
        control_var = control_var, snc_var = snc_var
      )

      # Same complete-case rows as the SNc-adjusted model, without the SNc term.
      same_dat <- adjusted_dat %>% select(-snc_raw, -snc_z)

      same_fits <- fit_prepared_model(same_dat)
      adj_fits  <- fit_prepared_model(adjusted_dat)

      if (inherits(same_fits$mixed, "error") || inherits(adj_fits$mixed, "error")) {
        return(tibble(
          predictor = predictor, outcome = outcome, snc_covariate = snc_covariate,
          status = "model error"
        ))
      }

      same <- extract_lmer_term(same_fits$mixed)
      adj  <- extract_lmer_term(adj_fits$mixed)
      
      same_lm <- extract_lm_term(same_fits$fixed)
      adj_lm  <- extract_lm_term(adj_fits$fixed)
      
      tibble(
        predictor = predictor,
        predictor_var = predictor_var,
        family = family,
        outcome = outcome,
        outcome_var = outcome_var,
        control_covariate = control_covariate,
        global_covariate = global_covariate,
        snc_covariate = snc_covariate,
        snc_qc_filter_applied = requires_snc_fw_qc,
        
        primary_n_hemi = primary_n_hemi,
        primary_n_sub = primary_n_sub,
        primary_beta = primary_beta,
        primary_p = primary_p,
        primary_q = primary_q,
        
        same_sample_n_hemi = nrow(same_dat),
        same_sample_n_sub = n_distinct(same_dat$sub_id),
        same_sample_beta = same$beta,
        same_sample_p = same$p,
        same_sample_lm_beta = same_lm$lm_beta,
        same_sample_lm_p = same_lm$lm_p,
        
        adjusted_n_hemi = nrow(adjusted_dat),
        adjusted_n_sub = n_distinct(adjusted_dat$sub_id),
        adjusted_beta = adj$beta,
        adjusted_p = adj$p,
        adjusted_lm_beta = adj_lm$lm_beta,
        adjusted_lm_p = adj_lm$lm_p,
        
        sample_restriction_beta_change = same$beta - primary_beta,
        snc_adjustment_beta_change = adj$beta - same$beta,
        total_beta_change = adj$beta - primary_beta,
        sign_reversal_after_snc_adjustment =
          sign(adj$beta) != sign(same$beta),
        
        singular = isSingular(adj_fits$mixed, tol = 1e-4),
        status = "OK"
      )
      
    })
  
} else {
  specificity_results <- tibble()
}



write_csv(
  specificity_results,
  file.path(out_dir, "section_2_3_snc_specificity_followups.csv")
)

cat("\nSNc regional-specificity follow-ups:\n")
print(specificity_results, n = Inf)

# ------------------------------------------------
# 9. Leave-one-participant-out sensitivity
# ------------------------------------------------

run_lopo <- function(effect_row) {

  full_dat <- prepare_model_data(
    long_all,
    effect_row$outcome_var,
    effect_row$predictor_var,
    effect_row$global_var,
    control_var = effect_row$control_var
  )

  participant_ids <- sort(unique(full_dat$sub_id))

  map_dfr(participant_ids, function(left_out_sub) {

    dat_i <- full_dat %>%
      filter(sub_id != left_out_sub)
      

    fits <- fit_prepared_model(dat_i)

    if (inherits(fits$mixed, "error")) {
      return(tibble(
        predictor = effect_row$predictor,
        outcome = effect_row$outcome,
        left_out_sub = left_out_sub,
        n_hemi = nrow(dat_i),
        n_sub = n_distinct(dat_i$sub_id),
        lopo_beta = NA_real_,
        lopo_p = NA_real_,
        sign_reversal = NA,
        p_lt_05 = NA,
        singular = NA,
        status = paste0("lmer error: ", fits$mixed$message)
      ))
    }

    stat <- extract_lmer_term(fits$mixed)

    tibble(
      predictor = effect_row$predictor,
      outcome = effect_row$outcome,
      left_out_sub = left_out_sub,
      n_hemi = nrow(dat_i),
      n_sub = n_distinct(dat_i$sub_id),
      original_beta = effect_row$primary_beta,
      original_p = effect_row$primary_p,
      original_q = effect_row$primary_q,
      lopo_beta = stat$beta,
      lopo_p = stat$p,
      sign_reversal = sign(stat$beta) != sign(effect_row$primary_beta),
      p_lt_05 = !is.na(stat$p) && stat$p < 0.05,
      singular = isSingular(fits$mixed, tol = 1e-4),
      status = "OK"
    )
  })
}

if (nrow(primary_effects) > 0) {
  lopo_results <- map_dfr(seq_len(nrow(primary_effects)), function(i) {
    run_lopo(primary_effects[i, ])
  })

  lopo_summary <- lopo_results %>%
    group_by(predictor, outcome, original_beta, original_p, original_q) %>%
    summarise(
      n_refits = n(),
      n_successful_refits = sum(status == "OK"),
      n_failed_refits = sum(status != "OK"),
      lopo_min_beta = min(lopo_beta, na.rm = TRUE),
      lopo_max_beta = max(lopo_beta, na.rm = TRUE),
      lopo_min_p = min(lopo_p, na.rm = TRUE),
      lopo_max_p = max(lopo_p, na.rm = TRUE),
      n_lopo_p_lt_05 = sum(p_lt_05, na.rm = TRUE),
      any_sign_reversal = any(sign_reversal, na.rm = TRUE),
      any_singular = any(singular, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(original_q, original_p)
} else {
  lopo_results <- tibble()
  lopo_summary <- tibble()
}

write_csv(lopo_results, file.path(out_dir, "section_2_3_leave_one_participant_results.csv"))
write_csv(lopo_summary, file.path(out_dir, "section_2_3_leave_one_participant_summary.csv"))

cat("\nLeave-one-participant-out summary:\n")
print(lopo_summary, n = Inf)

# ------------------------------------------------
# 10. Diagnostics and reproducibility information
# ------------------------------------------------

singular_primary <- ppn_primary_results %>%
  filter(!is.na(p), p < 0.05, singular) %>%
  select(predictor, outcome, beta, p, q_fdr, singular, lm_beta, lm_p)

write_csv(
  singular_primary,
  file.path(out_dir, "section_2_3_singular_primary_models_lm_comparison.csv")
)

complete_case_summary <- ppn_primary_results %>%
  group_by(predictor, family) %>%
  summarise(
    min_n_hemi = min(n_hemi, na.rm = TRUE),
    max_n_hemi = max(n_hemi, na.rm = TRUE),
    min_n_sub = min(n_sub, na.rm = TRUE),
    max_n_sub = max(n_sub, na.rm = TRUE),
    .groups = "drop"
  )

write_csv(
  complete_case_summary,
  file.path(out_dir, "section_2_3_complete_case_summary.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(out_dir, "sessionInfo.txt")
)

cat(
  "\nSection 2.3 PPN analysis complete.\n",
  "Outputs saved in: ", out_dir, "\n",
  sep = ""
)
