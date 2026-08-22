# Distinct nigral and brainstem pathology markers map onto separable subthalamic electrophysiological signatures in Parkinson's disease

This repository contains the de-identified derived data and R scripts used to reproduce the statistical analyses and principal figures reported in:

> Delgado-Sanchez A. et al. *Distinct nigral and brainstem pathology markers map onto separable subthalamic electrophysiological signatures in Parkinson's disease*.

The study combines OFF-medication subthalamic nucleus local field potential (STN LFP) recordings with quantitative MRI measures of substantia nigra pars compacta (SNc) and pedunculopontine nucleus (PPN) tissue properties in people with Parkinson's disease undergoing deep brain stimulation.

## Repository contents

```text
ADS_LFPvsMRI/
├── Analysis/
│   ├── 1.section_2_2_snc_primary.R
│   ├── 2.section_2_2_double_dissociation_bootstrap.R
│   ├── 3.section_2_3_ppn_primary.R
│   ├── 4.section_2_4_ppn_snc_moderation.R
│   └── 5.section_2_5_bradykinesia.R
├── data/
│   ├── combined_metrics.csv
│   └── qc_metrics.csv
├── Figures/
│   ├── Figure2.R
│   ├── Figure3.R
│   └── Figure4.R
└── README.md
```

The analysis scripts create a `results/` directory containing model results, sensitivity analyses and intermediate files. The figure scripts save publication-quality PDF and PNG files and the corresponding plot-level data.

## Data

`combined_metrics.csv` contains the de-identified derived variables required for the reported analyses. These include hemisphere-level STN electrophysiological features, regional and whole-brain MRI measures, relevant clinical measures and analysis covariates.

`qc_metrics.csv` contains hemisphere-level whole-SNc diffusion MRI voxel-retention metrics used for the SNc free-water quality-control procedure.

The shared analysis cohort comprises 35 participants. Individual analyses use model-specific complete-case samples, as described in the manuscript and recorded in the output files. Participant identifiers are pseudonymous. Users must not attempt to identify participants or link these data to other sources for that purpose.

## Software requirements

The analyses were written in R. Install the required packages with:

```r
install.packages(c(
  "lme4",
  "lmerTest",
  "dplyr",
  "tidyr",
  "purrr",
  "tibble",
  "readr",
  "car",
  "sandwich",
  "lmtest",
  "ggplot2",
  "patchwork",
  "stringr",
  "scales"
))
```

Package and R versions used for each analysis are written to `sessionInfo.txt` files in the relevant results directories.

## Running the analyses

Clone or download the repository and set the R working directory to its root:

```r
setwd("path/to/ADS_LFPvsMRI")
```

The scripts use relative paths and should be run in the following order:

```r
source("Analysis/1.section_2_2_snc_primary.R")
source("Analysis/2.section_2_2_double_dissociation_bootstrap.R")
source("Figures/Figure2.R")

source("Analysis/3.section_2_3_ppn_primary.R")
source("Analysis/4.section_2_4_ppn_snc_moderation.R")
source("Figures/Figure3.R")

source("Analysis/5.section_2_5_bradykinesia.R")
source("Figures/Figure4.R")
```

Figure 2 should be generated in the same R session immediately after Analysis Scripts 1 and 2 because it uses objects created by those scripts. Figures 3 and 4 read the saved analysis outputs and verify that plotted estimates agree with them.

The participant bootstrap in Script 2 uses a fixed random seed and 5,000 participant-level resamples. It may take appreciably longer than the other analyses.

## Analysis overview

### Section 2.2: nigral MRI markers and STN physiology

Script 1 tests associations between SNc free water or susceptibility and prespecified STN electrophysiological families. It applies the whole-SNc diffusion quality-control criterion, family-specific Benjamini-Hochberg correction, leave-one-participant-out sensitivity analyses and a stricter diffusion-QC sensitivity analysis.

Script 2 performs the matched-sample participant-level bootstrap test of the double-dissociation contrast between the two nigral MRI markers and aperiodic versus low-beta burst physiology.

### Section 2.3: PPN MRI markers and STN physiology

Script 3 tests PPN free-water-corrected axial diffusivity (cAD) and PPN free water as predictors of STN electrophysiology. Models include modality-matched whole-brain and ipsilateral corticospinal tract covariates. The script also performs regional-specificity and leave-one-participant-out analyses.

### Section 2.4: nigral moderation of PPN associations

Script 4 tests whether SNc susceptibility or free water moderates associations between PPN cAD and high-frequency aperiodic STN activity. Because the random-intercept models are singular, primary interaction inference uses linear models with participant-clustered HC2 robust standard errors.

### Section 2.5: electrophysiology and bradykinesia

Script 5 tests whether eight MRI-linked STN electrophysiological features are associated with contralateral OFF-medication bradykinesia. False-discovery-rate correction is applied separately to the four aperiodic and four low-beta features.

## Statistical conventions

- Analyses are performed at hemisphere level, with participant included as a random intercept where appropriate.
- Continuous variables are standardised within each model-specific complete-case sample.
- Leave-one-participant-out analyses remove both hemispheres together while retaining the standardisation derived from the original complete-case sample.
- All reported tests are two-sided.
- False discovery rate is controlled using the Benjamini-Hochberg procedure within the prespecified families described in the manuscript.

## Outputs

Running the scripts creates the following main directories:

```text
results/
├── section_2_2_snc_primary/
├── section_2_3_ppn/
├── section_2_4_ppn_snc_moderation/
├── section_2_5_bradykinesia/
└── figures/
```

Each analysis directory contains full model results, significant-result summaries, sensitivity analyses and a `sessionInfo.txt` file. Figure directories contain PDF and high-resolution PNG versions together with the data used for plotting.

## Reproducibility notes

Mixed-effects models with negligible participant-level random-intercept variance may produce `boundary (singular) fit` messages. These are expected and are explicitly assessed in the scripts. Corresponding fixed-effects results are saved where relevant. For the Section 2.4 interaction analyses, inference is based on participant-clustered HC2 robust standard errors, as specified in the manuscript.

## Citation

Please cite the associated manuscript when using this repository. A version-specific citation and DOI will be added when an archived release is created.

## Contact

Questions about the study or repository should be addressed to:

**Professor Nicola J. Ray**  
Manchester Metropolitan University  
[n.ray@mmu.ac.uk](mailto:n.ray@mmu.ac.uk)

## Funding

This work was supported by the Medical Research Council (MR/X005267/1).
