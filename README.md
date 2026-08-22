# ADS_LFPvsMRI
R scripts accompanying Delgado-Sanchez et al.: nigral and brainstem pathology associations with subthalamic electrophysiology

# Analysis code for Delgado-Sanchez et al.

This repository contains the R scripts used for the statistical analyses
reported in:

“Distinct nigral and brainstem pathology markers map onto separable
subthalamic electrophysiological signatures in Parkinson’s disease.”

## Included analyses

- `section_2_2_snc_primary.R`  
  Primary associations between substantia nigra pars compacta free water
  or susceptibility and subthalamic electrophysiological features.

- `section_2_2_double_dissociation_bootstrap.R`  
  Participant-level bootstrap analysis testing the double-dissociation
  contrast between nigral MRI markers and electrophysiological domains.

- `section_2_3_ppn.R`  
  Primary pedunculopontine nucleus analyses, corticospinal tract negative
  controls, regional-specificity analyses and sensitivity analyses.

- `section_2_4_ppn_snc_moderation.R`  
  Analyses testing whether nigral pathology moderates the association
  between pedunculopontine microstructure and high-frequency aperiodic
  subthalamic activity.

- `section_2_5_bradykinesia.R`  
  Analyses of associations between MRI-linked subthalamic
  electrophysiological features and contralateral OFF-medication
  bradykinesia.

## Data

The scripts use the following participant-level input files:

- `data/combined_metrics.csv`
- `data/qc_metrics.csv`

These data are not included in the repository. Access to
the study data is described in the Data Availability statement associated
with the manuscript.

The scripts are provided to document the analyses underlying the reported
findings. They call established publicly available R packages and are not
intended as a standalone software package.

## R packages

The scripts use established R packages including:

- readr
- dplyr
- tidyr
- purrr
- tibble
- stringr
- lme4
- lmerTest
- car
- sandwich
- lmtest

## Licence

The analysis scripts are released under the MIT License.
