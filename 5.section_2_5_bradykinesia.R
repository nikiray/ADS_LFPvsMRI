# ================================================================
# SECTION 2.5: MRI-linked STN electrophysiological features and
# contralateral OFF-medication bradykinesia
# ================================================================
# Eight electrophysiological features identified in Sections 2.2-2.4 are
# tested as predictors of contralateral OFF-medication bradykinesia:
# four aperiodic features; low-beta burst occupancy, median duration and rate;
# and low-beta oscillatory log spectral mass.
#
# Model:
#   z_bradykinesia ~ z_LFP + z_age + z_residual_OFF_LED + (1 | sub_id)
#
# Left STN features are paired with right-body bradykinesia and vice versa.
# Variables are standardised within each model-specific complete-case sample.
# BH FDR is applied separately across four aperiodic and four low-beta tests.
# The analysis cohort is restricted to participants contributing at least one
# MRI marker to the MRI-LFP study (32 of 35 participants). 
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
data_file <- file.path("data","combined_metrics.csv")
out_dir <- file.path("results","section_2_5_bradykinesia")


data_file <- "O:/Research/projects/MRC/Papers/shareableScripts/combined_metrics.csv"


dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)
led_covariate <- "rs-led_medoff"

# 2. Helpers
z_sample <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(x)
  s <- sd(x,na.rm=TRUE)
  if (is.na(s) || s==0) return(rep(NA_real_,length(x)))
  as.numeric((x-mean(x,na.rm=TRUE))/s)
}

check_cols <- function(dat,cols) {
  missing <- setdiff(cols,names(dat))
  if (length(missing)>0) stop("Missing expected columns:\n",paste(missing,collapse="\n"))
}

extract_lmer_term <- function(fit,term="predictor_z") {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) {
    return(tibble(beta=NA_real_,se=NA_real_,df=NA_real_,t=NA_real_,p=NA_real_))
  }
  tibble(
    beta=unname(sm[term,"Estimate"]),se=unname(sm[term,"Std. Error"]),
    df=if ("df" %in% colnames(sm)) unname(sm[term,"df"]) else NA_real_,
    t=unname(sm[term,"t value"]),p=unname(sm[term,"Pr(>|t|)"])
  )
}

extract_lm_term <- function(fit,term="predictor_z") {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) {
    return(tibble(lm_beta=NA_real_,lm_se=NA_real_,lm_t=NA_real_,lm_p=NA_real_))
  }
  tibble(
    lm_beta=unname(sm[term,"Estimate"]),lm_se=unname(sm[term,"Std. Error"]),
    lm_t=unname(sm[term,"t value"]),lm_p=unname(sm[term,"Pr(>|t|)"])
  )
}

get_max_vif <- function(fit) {
  v <- tryCatch(car::vif(fit),error=function(e) NA)
  if (length(v)==1 && all(is.na(v))) return(NA_real_)
  if (is.matrix(v) && "GVIF^(1/(2*Df))" %in% colnames(v)) {
    v <- v[,"GVIF^(1/(2*Df))"]
  } else if (is.matrix(v)) {
    v <- as.numeric(v)
  }
  v <- as.numeric(v)
  v <- v[is.finite(v)]
  if (length(v)==0) NA_real_ else max(v)
}

# 3. Read data and construct hemisphere-level dataset
wide <- read_csv(data_file,show_col_types=FALSE)

required_wide <- c(
  "sub_id","age",led_covariate,
  "right_updrs-bradykinesia_medoff","left_updrs-bradykinesia_medoff",
  "left_stn_4-39-fractal-offset_medsoff","right_stn_4-39-fractal-offset_medsoff",
  "left_stn_4-39-fractal-slope_medsoff","right_stn_4-39-fractal-slope_medsoff",
  "left_stn_40-80-fractal-offset_medsoff","right_stn_40-80-fractal-offset_medsoff",
  "left_stn_40-80-fractal-slope_medsoff","right_stn_40-80-fractal-slope_medsoff",
  "left_stn_13-20-thr-zrobust-occupancy_medoff",
  "right_stn_13-20-thr-zrobust-occupancy_medoff",
  "left_stn_13-20-thr-zrobust-bd-median_medoff",
  "right_stn_13-20-thr-zrobust-bd-median_medoff",
  "left_stn_13-20-thr-zrobust-burstrate_medoff",
  "right_stn_13-20-thr-zrobust-burstrate_medoff",
  "left_stn_13-20-oscillatory-logsm_medoff",
  "right_stn_13-20-oscillatory-logsm_medoff",
  "left_snc_fw","right_snc_fw","left_snc_chi","right_snc_chi",
  "left_ppn_fw","right_ppn_fw","left_ppn_fwcad","right_ppn_fwcad"
)
check_cols(wide,required_wide)

mri_study_cols <- c(
  "left_snc_fw","right_snc_fw","left_snc_chi","right_snc_chi",
  "left_ppn_fw","right_ppn_fw","left_ppn_fwcad","right_ppn_fwcad"
)

cohort_flags <- wide %>%
  transmute(
    sub_id,
    n_available_mri_markers=rowSums(!is.na(pick(all_of(mri_study_cols)))),
    include_mri_lfp_study_cohort=n_available_mri_markers>0
  )

analysis_ids <- cohort_flags %>%
  filter(include_mri_lfp_study_cohort) %>% pull(sub_id)
wide_analysis <- wide %>% filter(sub_id %in% analysis_ids)

write_csv(cohort_flags,file.path(out_dir,"section_2_5_participant_cohort_flags.csv"))

make_side <- function(dat,side=c("left","right")) {
  side <- match.arg(side)
  p <- paste0(side,"_")
  contralateral_brady_col <- if (side=="left") {
    "right_updrs-bradykinesia_medoff"
  } else {
    "left_updrs-bradykinesia_medoff"
  }
  dat %>% transmute(
    sub_id=sub_id,side=side,age=age,
    !!led_covariate := .data[[led_covariate]],
    brady_contralateral_off=.data[[contralateral_brady_col]],
    ap_4_39_offset=.data[[paste0(p,"stn_4-39-fractal-offset_medsoff")]],
    ap_4_39_slope=.data[[paste0(p,"stn_4-39-fractal-slope_medsoff")]],
    ap_40_80_offset=.data[[paste0(p,"stn_40-80-fractal-offset_medsoff")]],
    ap_40_80_slope=.data[[paste0(p,"stn_40-80-fractal-slope_medsoff")]],
    lb_occ=.data[[paste0(p,"stn_13-20-thr-zrobust-occupancy_medoff")]],
    lb_bd_median=.data[[paste0(p,"stn_13-20-thr-zrobust-bd-median_medoff")]],
    lb_burstrate=.data[[paste0(p,"stn_13-20-thr-zrobust-burstrate_medoff")]],
    low_beta_logsm=.data[[paste0(p,"stn_13-20-oscillatory-logsm_medoff")]]
  )
}

long_all <- bind_rows(
  make_side(wide_analysis,"left"),
  make_side(wide_analysis,"right")
)
lfp_available <- long_all %>%
  filter(if_any(c(ap_4_39_offset,ap_4_39_slope,ap_40_80_offset,
                  ap_40_80_slope,lb_occ,lb_bd_median,lb_burstrate,
                  low_beta_logsm),~!is.na(.x)))

cat(
  "\nParticipants in input cohort: ",n_distinct(wide$sub_id),
  "\nParticipants in MRI-LFP analysis cohort: ",n_distinct(wide_analysis$sub_id),
  "\nHemisphere rows constructed: ",nrow(long_all),
  "\nHemispheres with at least one carried-forward LFP feature: ",nrow(lfp_available),
  "\nParticipants with at least one carried-forward LFP feature: ",
  n_distinct(lfp_available$sub_id),"\n",sep=""
)
write_csv(long_all,file.path(out_dir,"section_2_5_hemisphere_level_data.csv"))

# 4. Prespecified LFP families
lfp_predictors <- tribble(
  ~fdr_family,~predictor_var,~predictor,~display_order,
  "aperiodic","ap_4_39_offset","4-39 Hz aperiodic offset",1,
  "aperiodic","ap_4_39_slope","4-39 Hz aperiodic slope",2,
  "aperiodic","ap_40_80_offset","40-80 Hz aperiodic offset",3,
  "aperiodic","ap_40_80_slope","40-80 Hz aperiodic slope",4,
  "low_beta","lb_occ","Low-beta percentage burst time",5,
  "low_beta","lb_bd_median","Low-beta median burst duration",6,
  "low_beta","lb_burstrate","Low-beta burst rate",7,
  "low_beta","low_beta_logsm","Low-beta oscillatory log spectral mass",8
)

# 5. Fit models
run_lfp_model <- function(fdr_family,predictor_var,predictor,display_order) {
  model_dat <- long_all %>%
    select(sub_id,side,age,all_of(led_covariate),
           brady_contralateral_off,all_of(predictor_var)) %>%
    rename(led_raw=all_of(led_covariate),
           outcome_raw=brady_contralateral_off,
           predictor_raw=all_of(predictor_var)) %>%
    drop_na(outcome_raw,predictor_raw,age,led_raw,sub_id) %>%
    mutate(outcome_z=z_sample(outcome_raw),predictor_z=z_sample(predictor_raw),
           age_z=z_sample(age),led_z=z_sample(led_raw))

  fit_mixed <- suppressWarnings(lmer(
    outcome_z ~ predictor_z + age_z + led_z + (1|sub_id),
    data=model_dat,REML=TRUE,
    control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5))
  ))
  fit_fixed <- lm(outcome_z ~ predictor_z + age_z + led_z,data=model_dat)
  mixed <- extract_lmer_term(fit_mixed)
  fixed <- extract_lm_term(fit_fixed)

  tibble(
    fdr_family=fdr_family,predictor=predictor,predictor_var=predictor_var,
    display_order=display_order,n_hemi=nrow(model_dat),
    n_sub=n_distinct(model_dat$sub_id),
    beta=mixed$beta,se=mixed$se,df=mixed$df,t=mixed$t,p=mixed$p,
    singular=isSingular(fit_mixed,tol=1e-4),
    lm_beta=fixed$lm_beta,lm_se=fixed$lm_se,lm_t=fixed$lm_t,lm_p=fixed$lm_p,
    max_vif=get_max_vif(fit_fixed)
  )
}

lfp_results <- pmap_dfr(lfp_predictors,run_lfp_model) %>%
  group_by(fdr_family) %>% mutate(q_fdr=p.adjust(p,method="BH")) %>%
  ungroup() %>%
  mutate(ci_lower_95=beta-1.96*se,ci_upper_95=beta+1.96*se,
         fdr_significant=!is.na(q_fdr) & q_fdr<0.05) %>%
  arrange(display_order)

write_csv(lfp_results,file.path(out_dir,"section_2_5_lfp_bradykinesia_all_models.csv"))
write_csv(lfp_results %>% filter(fdr_significant),
          file.path(out_dir,"section_2_5_lfp_bradykinesia_fdr_q_lt_0_05.csv"))
write_csv(lfp_results %>% filter(singular),
          file.path(out_dir,"section_2_5_singular_models_lm_comparison.csv"))

cat("\nSection 2.5 LFP-bradykinesia results:\n")
print(lfp_results %>%
        select(fdr_family,predictor,n_hemi,n_sub,beta,se,t,p,q_fdr,
               singular,lm_beta,lm_p,max_vif),n=Inf)

writeLines(capture.output(sessionInfo()),file.path(out_dir,"sessionInfo.txt"))
cat("\nSection 2.5 bradykinesia analysis complete.\n",
    "Outputs saved in: ",out_dir,"\n",sep="")
