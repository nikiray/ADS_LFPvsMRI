# ================================================================
# SECTION 2.2: SNc MRI predictors of OFF-medication STN physiology
# ================================================================
#
# Purpose
# -------
# Tests whether ipsilateral substantia nigra pars compacta (SNc) MRI markers
# are associated with OFF-medication STN electrophysiological features.
#
# Primary MRI predictors
# ----------------------
# SNc free water:
#   z_LFP ~ z_SNc_FW + z_age + z_residual_OFF_LED +
#           z_whole_brain_GM_FW + (1 | sub_id)
#
# SNc susceptibility:
#   z_LFP ~ z_SNc_susceptibility + z_age + z_residual_OFF_LED +
#           z_whole_brain_GM_susceptibility + (1 | sub_id)
#
# Key choices
# -----------
# - SNc free-water models use the whole-SNc DWI voxel-retention QC criterion.
# - The primary threshold requires >=10% of atlas-defined SNc voxels retained.
# - Continuous variables are z-standardised within each complete-case model.
# - FDR correction is applied separately within each MRI predictor and
#   electrophysiological family.
# - Corresponding fixed-effects models are retained to assess isolated
#   singular mixed-effects fits.
# - Leave-one-hemisphere-out sensitivity is run for nominally significant
#   MRI associations, retaining the original complete-case standardisation.
# ================================================================

library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(readr)
library(car)

# 1. User settings
data_file <- file.path("data", "combined_metrics.csv")
qc_file   <- file.path("data", "qc_metrics.csv")   # optional
out_dir <- file.path("results", "section_2_2_snc_primary")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

led_covariate <- "rs-led_medoff"
snc_qc_threshold <- 10
exclude_missing_snc_qc <- TRUE

# 2. Helper functions
z_sample <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(x)
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  as.numeric((x - mean(x, na.rm = TRUE)) / s)
}

check_cols <- function(dat, cols) {
  missing <- setdiff(cols, names(dat))
  if (length(missing) > 0) stop("Missing expected columns:\n", paste(missing, collapse = "\n"))
}

extract_lmer_term <- function(fit, term = "mri_z") {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) {
    return(tibble(beta=NA_real_, se=NA_real_, df=NA_real_, t=NA_real_, p=NA_real_))
  }
  tibble(
    beta = unname(sm[term, "Estimate"]),
    se   = unname(sm[term, "Std. Error"]),
    df   = if ("df" %in% colnames(sm)) unname(sm[term, "df"]) else NA_real_,
    t    = unname(sm[term, "t value"]),
    p    = unname(sm[term, "Pr(>|t|)"])
  )
}

extract_lm_term <- function(fit, term = "mri_z") {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) {
    return(tibble(lm_beta=NA_real_, lm_se=NA_real_, lm_t=NA_real_, lm_p=NA_real_))
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

# 3. Read data / construct hemisphere-level dataset
wide <- read_csv(data_file, show_col_types = FALSE)

if (file.exists(qc_file)) {
  qc <- read_csv(qc_file, show_col_types = FALSE)
  check_cols(qc, c("sub_id","left_snc_percvoxremaining","right_snc_percvoxremaining"))
  wide <- wide %>%
    select(-any_of(c("left_snc_percvoxremaining","right_snc_percvoxremaining"))) %>%
    left_join(
      qc %>% select(sub_id,left_snc_percvoxremaining,right_snc_percvoxremaining),
      by = "sub_id"
    )
}

required_wide <- c(
  "sub_id","age",led_covariate,
  "left_snc_fw","right_snc_fw",
  "left_snc_chi","right_snc_chi",
  "bilateral_cgm_fw","bilateral_cgm_chi",
  "left_snc_percvoxremaining","right_snc_percvoxremaining"
)
check_cols(wide, required_wide)

make_side <- function(dat, side=c("left","right")) {
  side <- match.arg(side)
  p <- paste0(side, "_")
  dat %>% transmute(
    sub_id=sub_id, side=side, age=age,
    !!led_covariate := .data[[led_covariate]],
    snc_fw=.data[[paste0(p,"snc_fw")]],
    snc_chi=.data[[paste0(p,"snc_chi")]],
    snc_percvoxremaining=.data[[paste0(p,"snc_percvoxremaining")]],
    bilateral_cgm_fw=bilateral_cgm_fw,
    bilateral_cgm_chi=bilateral_cgm_chi,
    ap_4_39_offset=.data[[paste0(p,"stn_4-39-fractal-offset_medsoff")]],
    ap_4_39_slope=.data[[paste0(p,"stn_4-39-fractal-slope_medsoff")]],
    ap_40_80_offset=.data[[paste0(p,"stn_40-80-fractal-offset_medsoff")]],
    ap_40_80_slope=.data[[paste0(p,"stn_40-80-fractal-slope_medsoff")]],
    osc_13_20_logsm=.data[[paste0(p,"stn_13-20-oscillatory-logsm_medoff")]],
    osc_21_30_logsm=.data[[paste0(p,"stn_21-30-oscillatory-logsm_medoff")]],
    lb_occ=.data[[paste0(p,"stn_13-20-thr-zrobust-occupancy_medoff")]],
    lb_burstrate=.data[[paste0(p,"stn_13-20-thr-zrobust-burstrate_medoff")]],
    lb_bd_median=.data[[paste0(p,"stn_13-20-thr-zrobust-bd-median_medoff")]],
    lb_ibi_median=.data[[paste0(p,"stn_13-20-thr-zrobust-ibi-median_medoff")]],
    hb_occ=.data[[paste0(p,"stn_21-30-thr-zrobust-occupancy_medoff")]],
    hb_burstrate=.data[[paste0(p,"stn_21-30-thr-zrobust-burstrate_medoff")]],
    hb_bd_median=.data[[paste0(p,"stn_21-30-thr-zrobust-bd-median_medoff")]],
    hb_ibi_median=.data[[paste0(p,"stn_21-30-thr-zrobust-ibi-median_medoff")]]
  )
}

long_all <- bind_rows(make_side(wide,"left"), make_side(wide,"right")) %>%
  mutate(
    snc_fw_qc_pass = case_when(
      is.na(snc_percvoxremaining) & exclude_missing_snc_qc ~ FALSE,
      is.na(snc_percvoxremaining) & !exclude_missing_snc_qc ~ TRUE,
      snc_percvoxremaining >= snc_qc_threshold ~ TRUE,
      TRUE ~ FALSE
    )
  )

# Export the pre-QC hemisphere-level dataset required by the
# Section 2.2 double-dissociation bootstrap script.
write_csv(
  long_all,
  file.path(out_dir, "section_2_2_hemisphere_long_all_preQC.csv")
)

long_fw_qc <- long_all %>% filter(snc_fw_qc_pass)

qc_exclusions <- long_all %>%
  filter(!snc_fw_qc_pass) %>%
  select(sub_id, side, snc_percvoxremaining)

write_csv(qc_exclusions, file.path(out_dir,"section_2_2_snc_fw_qc_exclusions.csv"))

cat(
  "\nHemisphere rows before SNc DWI QC:", nrow(long_all),
  "\nParticipants before SNc DWI QC:", n_distinct(long_all$sub_id),
  "\nHemisphere rows after >=10% SNc voxel-retention QC:", nrow(long_fw_qc),
  "\nParticipants after SNc DWI QC:", n_distinct(long_fw_qc$sub_id),
  "\n"
)

# 4. Predictors / outcome families
mri_predictors <- tribble(
  ~predictor_var, ~predictor, ~global_var, ~global_covariate, ~requires_dwi_qc,
  "snc_fw","SNc free water","bilateral_cgm_fw","whole-brain GM free water",TRUE,
  "snc_chi","SNc susceptibility","bilateral_cgm_chi","whole-brain GM susceptibility",FALSE
)

lfp_outcomes <- tribble(
  ~outcome_var, ~outcome, ~family,
  "ap_4_39_offset","4-39 Hz aperiodic offset","aperiodic",
  "ap_4_39_slope","4-39 Hz aperiodic slope","aperiodic",
  "ap_40_80_offset","40-80 Hz aperiodic offset","aperiodic",
  "ap_40_80_slope","40-80 Hz aperiodic slope","aperiodic",
  "osc_13_20_logsm","Low-beta oscillatory log spectral mass","low_beta",
  "lb_occ","Low-beta percentage burst time","low_beta",
  "lb_burstrate","Low-beta burst rate","low_beta",
  "lb_bd_median","Low-beta median burst duration","low_beta",
  "lb_ibi_median","Low-beta median inter-burst interval","low_beta",
  "osc_21_30_logsm","High-beta oscillatory log spectral mass","high_beta",
  "hb_occ","High-beta percentage burst time","high_beta",
  "hb_burstrate","High-beta burst rate","high_beta",
  "hb_bd_median","High-beta median burst duration","high_beta",
  "hb_ibi_median","High-beta median inter-burst interval","high_beta"
)

# 5. Model functions
prepare_model_data <- function(dat,outcome_var,predictor_var,global_var) {
  cols <- c("sub_id","side","age",led_covariate,outcome_var,predictor_var,global_var)
  check_cols(dat, cols)
  dat %>%
    select(all_of(cols)) %>%
    rename(
      outcome_raw=all_of(outcome_var),
      mri_raw=all_of(predictor_var),
      global_raw=all_of(global_var),
      led_raw=all_of(led_covariate)
    ) %>%
    drop_na(outcome_raw,mri_raw,global_raw,age,led_raw,sub_id) %>%
    mutate(
      outcome_z=z_sample(outcome_raw),
      mri_z=z_sample(mri_raw),
      global_z=z_sample(global_raw),
      age_z=z_sample(age),
      led_z=z_sample(led_raw)
    )
}

fit_models <- function(dat) {
  f_mixed <- outcome_z ~ mri_z + age_z + led_z + global_z + (1|sub_id)
  f_fixed <- outcome_z ~ mri_z + age_z + led_z + global_z
  list(
    mixed=suppressWarnings(lmer(
      f_mixed,data=dat,REML=TRUE,
      control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
    )),
    fixed=lm(f_fixed,data=dat)
  )
}

run_leave_one_hemisphere <- function(model_dat, original_beta) {
  map_dfr(seq_len(nrow(model_dat)), function(i) {
    dat_i <- model_dat[-i,,drop=FALSE]  # retain original scaling
    fit_i <- tryCatch(
      suppressWarnings(lmer(
        outcome_z ~ mri_z + age_z + led_z + global_z + (1|sub_id),
        data=dat_i,REML=TRUE,
        control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
      )),
      error=function(e) e
    )
    if (inherits(fit_i,"error")) {
      return(tibble(
        left_out_sub=model_dat$sub_id[i],
        left_out_side=model_dat$side[i],
        loo_beta=NA_real_, loo_p=NA_real_,
        sign_reversal=NA, status=paste0("lmer error: ",fit_i$message)
      ))
    }
    stat <- extract_lmer_term(fit_i)
    tibble(
      left_out_sub=model_dat$sub_id[i],
      left_out_side=model_dat$side[i],
      loo_beta=stat$beta,
      loo_p=stat$p,
      sign_reversal=sign(stat$beta)!=sign(original_beta),
      status="OK"
    )
  })
}

run_one_model <- function(dat,outcome_var,outcome,family,predictor_var,predictor,
                          global_var,global_covariate,dwi_qc_filter_applied) {
  model_dat <- prepare_model_data(dat,outcome_var,predictor_var,global_var)
  fits <- fit_models(model_dat)
  mixed <- extract_lmer_term(fits$mixed)
  fixed <- extract_lm_term(fits$fixed)

  out <- tibble(
    predictor=predictor, predictor_var=predictor_var,
    outcome=outcome, outcome_var=outcome_var, family=family,
    dwi_qc_filter_applied=dwi_qc_filter_applied,
    global_covariate=global_covariate,
    n_hemi=nrow(model_dat), n_sub=n_distinct(model_dat$sub_id),
    beta=mixed$beta, se=mixed$se, df=mixed$df, t=mixed$t, p=mixed$p,
    singular=isSingular(fits$mixed,tol=1e-4),
    lm_beta=fixed$lm_beta, lm_se=fixed$lm_se, lm_t=fixed$lm_t, lm_p=fixed$lm_p,
    max_vif=get_max_vif(fits$fixed)
  )

  if (!is.na(mixed$p) && mixed$p < 0.05) {
    loo <- run_leave_one_hemisphere(model_dat,mixed$beta)
    write_csv(
      loo,
      file.path(out_dir,paste0("section_2_2_",predictor_var,"__",outcome_var,"_leave_one_hemisphere.csv"))
    )
    out <- out %>% mutate(
      loo_min_beta=min(loo$loo_beta,na.rm=TRUE),
      loo_max_beta=max(loo$loo_beta,na.rm=TRUE),
      loo_max_p=max(loo$loo_p,na.rm=TRUE),
      loo_sign_reversal=any(loo$sign_reversal,na.rm=TRUE)
    )
  } else {
    out <- out %>% mutate(
      loo_min_beta=NA_real_, loo_max_beta=NA_real_,
      loo_max_p=NA_real_, loo_sign_reversal=NA
    )
  }
  out
}

# 6. Fit all primary models / FDR
results <- crossing(mri_predictors,lfp_outcomes) %>%
  pmap_dfr(function(predictor_var,predictor,global_var,global_covariate,requires_dwi_qc,
                    outcome_var,outcome,family) {
    analysis_dat <- if (requires_dwi_qc) long_fw_qc else long_all
    run_one_model(
      analysis_dat,outcome_var,outcome,family,predictor_var,predictor,
      global_var,global_covariate,requires_dwi_qc
    )
  }) %>%
  group_by(predictor,family) %>%
  mutate(q_fdr=p.adjust(p,method="BH")) %>%
  ungroup() %>%
  arrange(predictor,family,q_fdr,p)

fdr_significant <- results %>%
  filter(!is.na(q_fdr),q_fdr < 0.05) %>%
  arrange(q_fdr,p)

write_csv(results,file.path(out_dir,"section_2_2_snc_primary_all_models.csv"))
write_csv(fdr_significant,file.path(out_dir,"section_2_2_snc_primary_fdr_q_lt_0_05.csv"))

cat("\nFDR-significant Section 2.2 SNc results:\n")
print(
  fdr_significant %>%
    select(
      predictor,outcome,n_hemi,n_sub,beta,se,t,p,q_fdr,
      singular,lm_beta,lm_p,loo_min_beta,loo_max_beta,
      loo_max_p,loo_sign_reversal
    ),
  n=Inf
)

# 7. Checks / reproducibility
singular_nominal <- results %>%
  filter(!is.na(p),p < 0.05,singular) %>%
  select(predictor,outcome,beta,p,q_fdr,lm_beta,lm_p)

write_csv(
  singular_nominal,
  file.path(out_dir,"section_2_2_singular_nominal_models_lm_comparison.csv")
)

complete_case_summary <- results %>%
  group_by(predictor,family,dwi_qc_filter_applied) %>%
  summarise(
    min_n_hemi=min(n_hemi,na.rm=TRUE),
    max_n_hemi=max(n_hemi,na.rm=TRUE),
    min_n_sub=min(n_sub,na.rm=TRUE),
    max_n_sub=max(n_sub,na.rm=TRUE),
    .groups="drop"
  )

write_csv(
  complete_case_summary,
  file.path(out_dir,"section_2_2_complete_case_summary.csv")
)

writeLines(capture.output(sessionInfo()), file.path(out_dir,"sessionInfo.txt"))

cat("\nSection 2.2 SNc primary analysis complete.\n",
    "Outputs saved in: ", out_dir, "\n", sep="")
