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
# - FDR-significant associations undergo leave-one-participant-out (LOPO)
#   sensitivity analysis. Both hemispheres are removed together and the
#   standardisation from the original complete-case model is retained.
# - FDR-significant SNc free-water associations are repeated using the
#   prespecified stricter >=30% voxel-retention threshold.
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
snc_qc_threshold_strict <- 30
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

run_leave_one_participant <- function(model_dat, original_beta) {
  participant_ids <- sort(unique(model_dat$sub_id))

  map_dfr(participant_ids, function(left_out_sub) {
    # model_dat was standardised once using the original complete-case sample.
    # Filtering it here therefore retains the original standardisation.
    dat_i <- model_dat %>% filter(sub_id != left_out_sub)
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
        left_out_sub=left_out_sub,
        n_hemi=nrow(dat_i), n_sub=n_distinct(dat_i$sub_id),
        lopo_beta=NA_real_, lopo_p=NA_real_,
        sign_reversal=NA, status=paste0("lmer error: ",fit_i$message)
      ))
    }
    stat <- extract_lmer_term(fit_i)
    tibble(
      left_out_sub=left_out_sub,
      n_hemi=nrow(dat_i), n_sub=n_distinct(dat_i$sub_id),
      lopo_beta=stat$beta,
      lopo_p=stat$p,
      sign_reversal=sign(stat$beta)!=sign(original_beta),
      p_lt_05=!is.na(stat$p) && stat$p < 0.05,
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
    global_covariate=global_covariate, global_var=global_var,
    n_hemi=nrow(model_dat), n_sub=n_distinct(model_dat$sub_id),
    beta=mixed$beta, se=mixed$se, df=mixed$df, t=mixed$t, p=mixed$p,
    singular=isSingular(fits$mixed,tol=1e-4),
    lm_beta=fixed$lm_beta, lm_se=fixed$lm_se, lm_t=fixed$lm_t, lm_p=fixed$lm_p,
    max_vif=get_max_vif(fits$fixed)
  )

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
      singular,lm_beta,lm_p,max_vif
    ),
  n=Inf
)

# 7. Leave-one-participant-out sensitivity for FDR-significant findings
lopo_results <- map_dfr(seq_len(nrow(fdr_significant)), function(i) {
  rr <- fdr_significant[i, ]
  analysis_dat <- if (rr$dwi_qc_filter_applied[[1]]) long_fw_qc else long_all
  model_dat <- prepare_model_data(
    analysis_dat,
    rr$outcome_var[[1]],
    rr$predictor_var[[1]],
    rr$global_var[[1]]
  )

  run_leave_one_participant(model_dat, rr$beta[[1]]) %>%
    mutate(
      predictor=rr$predictor[[1]], outcome=rr$outcome[[1]],
      original_beta=rr$beta[[1]], original_p=rr$p[[1]],
      original_q=rr$q_fdr[[1]], .before=1
    )
})

lopo_summary <- lopo_results %>%
  group_by(predictor,outcome,original_beta,original_p,original_q) %>%
  summarise(
    n_refits=n(),
    n_successful_refits=sum(status=="OK"),
    n_failed_refits=sum(status!="OK"),
    lopo_min_beta=ifelse(all(is.na(lopo_beta)),NA_real_,min(lopo_beta,na.rm=TRUE)),
    lopo_max_beta=ifelse(all(is.na(lopo_beta)),NA_real_,max(lopo_beta,na.rm=TRUE)),
    n_lopo_p_lt_05=sum(p_lt_05,na.rm=TRUE),
    any_sign_reversal=any(sign_reversal,na.rm=TRUE),
    .groups="drop"
  ) %>% arrange(original_q,original_p)

write_csv(lopo_results,file.path(out_dir,"section_2_2_leave_one_participant_results.csv"))
write_csv(lopo_summary,file.path(out_dir,"section_2_2_leave_one_participant_summary.csv"))

cat("\nLeave-one-participant-out summary for FDR-significant findings:\n")
print(lopo_summary,n=Inf)

# 8. Stricter 30% SNc voxel-retention sensitivity
long_fw_qc_strict <- long_all %>%
  filter(!is.na(snc_percvoxremaining),
         snc_percvoxremaining >= snc_qc_threshold_strict)

strict_fw_effects <- fdr_significant %>%
  filter(predictor_var=="snc_fw")

strict_qc_results <- map_dfr(seq_len(nrow(strict_fw_effects)), function(i) {
  rr <- strict_fw_effects[i, ]
  run_one_model(
    long_fw_qc_strict,
    rr$outcome_var[[1]],rr$outcome[[1]],rr$family[[1]],
    rr$predictor_var[[1]],rr$predictor[[1]],
    rr$global_var[[1]],rr$global_covariate[[1]],TRUE
  ) %>%
    mutate(
      qc_threshold_percent=snc_qc_threshold_strict,
      primary_beta=rr$beta[[1]], primary_p=rr$p[[1]],
      primary_q=rr$q_fdr[[1]], .before=1
    )
})

write_csv(strict_qc_results,file.path(out_dir,"section_2_2_snc_fw_strict_30pct_qc_sensitivity.csv"))
cat("\nStricter >=30% SNc voxel-retention sensitivity:\n")
print(strict_qc_results %>% select(predictor,outcome,n_hemi,n_sub,beta,se,t,p),n=Inf)

# 9. Checks / reproducibility
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
