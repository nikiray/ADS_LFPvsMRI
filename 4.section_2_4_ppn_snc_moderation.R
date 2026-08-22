# ================================================================
# SECTION 2.4: Nigral moderation of PPN cAD associations with
# high-frequency aperiodic STN activity

# ================================================================
#
# Purpose
# -------
# Reproduces the Section 2.4 moderation analyses:
#   1. PPN cAD x SNc susceptibility interactions
#   2. PPN cAD x SNc free-water interactions
#   3. CST cAD x SNc-marker negative-control interactions
#   4. Leave-one-participant-out (LOPO) sensitivity analyses
#
# Primary interaction model:
#   z_LFP ~ z_PPN_cAD * z_SNc_marker +
#           z_CST_cAD +
#           z_age + z_residual_OFF_LED + z_whole_brain_cAD
#
# Participant-clustered HC2 standard errors from the corresponding lm model
# are used for the primary interaction inference because the mixed-effects
# models are expected to be singular when participant-level random-intercept
# variance is negligible.
#
# Key choices
# -----------
# - Outcomes: 40-80 Hz aperiodic offset and slope.
# - Moderators: ipsilateral SNc susceptibility and SNc free water.
# - The 10% voxel-retention threshold is applied
# - BH correction is applied across the four primary PPN interaction tests.
# - CST negative-control interactions are corrected separately.
# - LOPO removes both hemispheres for one participant at a time while retaining
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

if (!requireNamespace("sandwich", quietly = TRUE) ||
    !requireNamespace("lmtest", quietly = TRUE)) {
  stop("Packages 'sandwich' and 'lmtest' are required.")
}

# 1. User settings
data_file <- file.path("data", "combined_metrics.csv")
qc_file   <- file.path("data", "qc_metrics.csv")  # optional


out_dir <- file.path("results", "section_2_4_ppn_snc_moderation")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
led_covariate <- "rs-led_medoff"
snc_qc_threshold <- 10
exclude_missing_snc_qc <- TRUE

# 2. Helpers
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

extract_lmer_term <- function(fit, term) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(tibble(beta=NA_real_, se=NA_real_, df=NA_real_, t=NA_real_, p=NA_real_))
  tibble(
    beta = unname(sm[term, "Estimate"]),
    se   = unname(sm[term, "Std. Error"]),
    df   = if ("df" %in% colnames(sm)) unname(sm[term, "df"]) else NA_real_,
    t    = unname(sm[term, "t value"]),
    p    = unname(sm[term, "Pr(>|t|)"])
  )
}

extract_lm_term <- function(fit, term) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) return(tibble(lm_beta=NA_real_, lm_se=NA_real_, lm_t=NA_real_, lm_p=NA_real_))
  tibble(
    lm_beta = unname(sm[term, "Estimate"]),
    lm_se   = unname(sm[term, "Std. Error"]),
    lm_t    = unname(sm[term, "t value"]),
    lm_p    = unname(sm[term, "Pr(>|t|)"])
  )
}

extract_clustered_term <- function(fit_lm, model_dat, term="predictor_z:moderator_z") {
  vc <- sandwich::vcovCL(fit_lm, cluster=model_dat$sub_id, type="HC2")
  ct <- lmtest::coeftest(fit_lm, vcov.=vc)
  if (!term %in% rownames(ct)) {
    return(tibble(cluster_beta=NA_real_, cluster_se=NA_real_, cluster_t=NA_real_,
                  cluster_p_two_tailed=NA_real_, cluster_status="Interaction term missing"))
  }
  tibble(
    cluster_beta = unname(ct[term,1]),
    cluster_se = unname(ct[term,2]),
    cluster_t = unname(ct[term,3]),
    cluster_p_two_tailed = unname(ct[term,4]),
    cluster_status = "OK"
  )
}

safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm=TRUE)
safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm=TRUE)

# 3. Read data / hemisphere dataset
wide <- read_csv(data_file, show_col_types = FALSE)

if (file.exists(qc_file)) {
  qc <- read_csv(qc_file, show_col_types = FALSE)
  check_cols(qc, c("sub_id","left_snc_percvoxremaining","right_snc_percvoxremaining"))
  wide <- wide %>%
    select(-any_of(c("left_snc_percvoxremaining","right_snc_percvoxremaining"))) %>%
    left_join(qc %>% select(sub_id,left_snc_percvoxremaining,right_snc_percvoxremaining), by="sub_id")
}

required_wide <- c(
  "sub_id","age",led_covariate,
  "left_ppn_fwcad","right_ppn_fwcad",
  "left_cst_fwcad","right_cst_fwcad",
  "bilateral_cgm_fwcad",
  "left_snc_fw","right_snc_fw",
  "left_snc_chi","right_snc_chi",
  "left_snc_percvoxremaining","right_snc_percvoxremaining"
)
check_cols(wide, required_wide)

make_side <- function(dat, side=c("left","right")) {
  side <- match.arg(side)
  p <- paste0(side, "_")
  dat %>% transmute(
    sub_id=sub_id, side=side, age=age,
    !!led_covariate := .data[[led_covariate]],
    ppn_ad=.data[[paste0(p,"ppn_fwcad")]],
    cst_ad=.data[[paste0(p,"cst_fwcad")]],
    bilateral_cgm_ad=bilateral_cgm_fwcad,
    snc_fw=.data[[paste0(p,"snc_fw")]],
    snc_chi=.data[[paste0(p,"snc_chi")]],
    snc_percvoxremaining=.data[[paste0(p,"snc_percvoxremaining")]],
    ap_40_80_offset=.data[[paste0(p,"stn_40-80-fractal-offset_medsoff")]],
    ap_40_80_slope=.data[[paste0(p,"stn_40-80-fractal-slope_medsoff")]]
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

# 4. Specifications
interaction_outcomes <- tribble(
  ~outcome_var, ~outcome,
  "ap_40_80_offset","40-80 Hz aperiodic offset",
  "ap_40_80_slope","40-80 Hz aperiodic slope"
)

primary_specs <- tribble(
  ~analysis_role, ~interaction, ~predictor_var, ~predictor,
  ~moderator_var, ~moderator, ~control_var, ~control_covariate,
  ~global_var, ~global_covariate, ~apply_snc_dwi_qc,
  "primary","PPN cAD x SNc susceptibility","ppn_ad","PPN cAD","snc_chi","SNc susceptibility",
  "cst_ad","ipsilateral CST cAD","bilateral_cgm_ad","whole-brain GM cAD",FALSE,
  "primary","PPN cAD x SNc free water","ppn_ad","PPN cAD","snc_fw","SNc free water",
  "cst_ad","ipsilateral CST cAD","bilateral_cgm_ad","whole-brain GM cAD",TRUE
)

negative_control_specs <- tribble(
  ~analysis_role, ~interaction, ~predictor_var, ~predictor,
  ~moderator_var, ~moderator, ~control_var, ~control_covariate,
  ~global_var, ~global_covariate, ~apply_snc_dwi_qc,
  "negative_control","CST cAD x SNc susceptibility","cst_ad","CST cAD","snc_chi","SNc susceptibility",
  NA_character_,NA_character_,"bilateral_cgm_ad","whole-brain GM cAD",FALSE,
  "negative_control","CST cAD x SNc free water","cst_ad","CST cAD","snc_fw","SNc free water",
  NA_character_,NA_character_,"bilateral_cgm_ad","whole-brain GM cAD",TRUE
)

all_specs <- bind_rows(primary_specs, negative_control_specs)

# 5. Model functions
make_interaction_data <- function(dat,outcome_var,predictor_var,moderator_var,global_var,control_var=NULL) {
  cols <- c("sub_id","side","age",led_covariate,outcome_var,predictor_var,moderator_var,global_var,control_var)
  cols <- cols[!is.na(cols) & nzchar(cols)]
  out <- dat %>% select(all_of(cols)) %>%
    rename(outcome_raw=all_of(outcome_var), predictor_raw=all_of(predictor_var),
           moderator_raw=all_of(moderator_var), global_raw=all_of(global_var),
           led_raw=all_of(led_covariate))
  if (!is.null(control_var)) out <- out %>% rename(control_raw=all_of(control_var))
  req <- c("outcome_raw","predictor_raw","moderator_raw","global_raw","age","led_raw","sub_id")
  if (!is.null(control_var)) req <- c(req,"control_raw")
  out <- out %>% drop_na(all_of(req)) %>%
    mutate(
      outcome_z=z_sample(outcome_raw), predictor_z=z_sample(predictor_raw),
      moderator_z=z_sample(moderator_raw), global_z=z_sample(global_raw),
      age_z=z_sample(age), led_z=z_sample(led_raw)
    )
  if (!is.null(control_var)) out <- out %>% mutate(control_z=z_sample(control_raw))
  out
}

fit_interaction_models <- function(dat) {
  if ("control_z" %in% names(dat)) {
    f_lmer <- outcome_z ~ predictor_z*moderator_z + control_z + age_z + led_z + global_z + (1|sub_id)
    f_lm   <- outcome_z ~ predictor_z*moderator_z + control_z + age_z + led_z + global_z
  } else {
    f_lmer <- outcome_z ~ predictor_z*moderator_z + age_z + led_z + global_z + (1|sub_id)
    f_lm   <- outcome_z ~ predictor_z*moderator_z + age_z + led_z + global_z
  }
  list(
    mixed = suppressWarnings(lmer(f_lmer,data=dat,REML=TRUE,
                                  control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))),
    fixed = lm(f_lm,data=dat)
  )
}

run_interaction_model <- function(dat,analysis_role,interaction,predictor_var,predictor,
                                  moderator_var,moderator,control_var,control_covariate,
                                  global_var,global_covariate,apply_snc_dwi_qc,outcome_var,outcome) {
  analysis_dat <- if (apply_snc_dwi_qc) dat %>% filter(snc_fw_qc_pass) else dat
  control_arg <- if (is.na(control_var)) NULL else control_var
  model_dat <- make_interaction_data(analysis_dat,outcome_var,predictor_var,moderator_var,global_var,control_arg)
  fits <- fit_interaction_models(model_dat)
  term <- "predictor_z:moderator_z"
  mixed <- extract_lmer_term(fits$mixed,term)
  lm_stats <- extract_lm_term(fits$fixed,term)
  clustered <- extract_clustered_term(fits$fixed,model_dat,term)
  tibble(
    analysis_role=analysis_role, interaction=interaction,
    predictor=predictor, predictor_var=predictor_var,
    moderator=moderator, moderator_var=moderator_var,
    control_covariate=control_covariate, control_var=control_var,
    global_covariate=global_covariate, global_var=global_var,
    outcome=outcome, outcome_var=outcome_var,
    snc_dwi_qc_applied=apply_snc_dwi_qc,
    n_hemi=nrow(model_dat), n_sub=n_distinct(model_dat$sub_id),
    mixed_beta=mixed$beta, mixed_se=mixed$se, mixed_t=mixed$t,
    mixed_p_two_tailed=mixed$p,
    singular=isSingular(fits$mixed,tol=1e-4),
    lm_beta=lm_stats$lm_beta, lm_p_two_tailed=lm_stats$lm_p
  ) %>% bind_cols(clustered)
}

# 6. Run models / FDR
interaction_results <- crossing(all_specs, interaction_outcomes) %>%
  pmap_dfr(function(analysis_role,interaction,predictor_var,predictor,moderator_var,moderator,
                    control_var,control_covariate,global_var,global_covariate,apply_snc_dwi_qc,
                    outcome_var,outcome) {
    run_interaction_model(long_all,analysis_role,interaction,predictor_var,predictor,
                          moderator_var,moderator,control_var,control_covariate,
                          global_var,global_covariate,apply_snc_dwi_qc,outcome_var,outcome)
  })

primary_results <- interaction_results %>%
  filter(analysis_role=="primary") %>%
  mutate(cluster_q_primary_family=p.adjust(cluster_p_two_tailed,method="BH"))

negative_control_results <- interaction_results %>%
  filter(analysis_role=="negative_control") %>%
  mutate(cluster_q_negative_control=p.adjust(cluster_p_two_tailed,method="BH"))

interaction_results <- interaction_results %>%
  left_join(primary_results %>% select(interaction,outcome,cluster_q_primary_family),
            by=c("interaction","outcome")) %>%
  left_join(negative_control_results %>% select(interaction,outcome,cluster_q_negative_control),
            by=c("interaction","outcome"))

write_csv(interaction_results,file.path(out_dir,"section_2_4_interaction_all_models.csv"))
write_csv(interaction_results %>% filter(analysis_role=="primary"),
          file.path(out_dir,"section_2_4_primary_ppn_interactions.csv"))
write_csv(interaction_results %>% filter(analysis_role=="negative_control"),
          file.path(out_dir,"section_2_4_cst_negative_control_interactions.csv"))

cat("\nPrimary PPN interaction results:\n")
print(interaction_results %>% filter(analysis_role=="primary") %>%
        select(interaction,outcome,n_hemi,n_sub,cluster_beta,cluster_se,cluster_t,
               cluster_p_two_tailed,cluster_q_primary_family,mixed_beta,
               mixed_p_two_tailed,singular,lm_beta,lm_p_two_tailed), n=Inf)

# 7. LOPO
run_lopo_interaction <- function(rr) {
  analysis_dat <- if (rr$snc_dwi_qc_applied[[1]]) long_all %>% filter(snc_fw_qc_pass) else long_all
  control_arg <- if (is.na(rr$control_var[[1]])) NULL else rr$control_var[[1]]
  full_dat <- make_interaction_data(
    analysis_dat, rr$outcome_var[[1]], rr$predictor_var[[1]],
    rr$moderator_var[[1]], rr$global_var[[1]], control_arg
  )
  ids <- sort(unique(full_dat$sub_id))
  map_dfr(ids,function(id){
    dat_i <- full_dat %>% filter(sub_id != id)  # retain original scaling
    fits <- fit_interaction_models(dat_i)
    cl <- extract_clustered_term(fits$fixed,dat_i)
    tibble(
      interaction=rr$interaction[[1]], outcome=rr$outcome[[1]],
      left_out_sub=id, n_hemi=nrow(dat_i), n_sub=n_distinct(dat_i$sub_id),
      original_cluster_beta=rr$cluster_beta[[1]],
      lopo_cluster_beta=cl$cluster_beta,
      lopo_cluster_p=cl$cluster_p_two_tailed,
      sign_reversal=sign(cl$cluster_beta)!=sign(rr$cluster_beta[[1]]),
      p_lt_05=cl$cluster_p_two_tailed < .05,
      status=cl$cluster_status
    )
  })
}

primary_full <- interaction_results %>% filter(analysis_role=="primary")
lopo_results <- map_dfr(seq_len(nrow(primary_full)), function(i) run_lopo_interaction(primary_full[i,]))

lopo_summary <- lopo_results %>%
  group_by(interaction,outcome,original_cluster_beta) %>%
  summarise(
    n_refits=n(),
    n_successful_refits=sum(status=="OK"),
    n_failed_refits=sum(status!="OK"),
    lopo_min_beta=safe_min(lopo_cluster_beta),
    lopo_max_beta=safe_max(lopo_cluster_beta),
    lopo_min_p=safe_min(lopo_cluster_p),
    lopo_max_p=safe_max(lopo_cluster_p),
    n_lopo_p_lt_05=sum(p_lt_05,na.rm=TRUE),
    any_sign_reversal=any(sign_reversal,na.rm=TRUE),
    .groups="drop"
  )

write_csv(lopo_results,file.path(out_dir,"section_2_4_leave_one_participant_interactions.csv"))
write_csv(lopo_summary,file.path(out_dir,"section_2_4_leave_one_participant_interaction_summary.csv"))

cat("\nLOPO summary:\n")
print(lopo_summary,n=Inf)

# 8. Reproducibility
writeLines(capture.output(sessionInfo()), file.path(out_dir,"sessionInfo.txt"))

cat("\nSection 2.4 moderation analysis complete.\n",
    "Outputs saved in: ", out_dir, "\n", sep="")
