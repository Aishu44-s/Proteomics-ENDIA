# ENDIA Proteomics Analysis

<p align="center">
  <strong>Reproducible R workflows for maternal proteomics analyses in the ENDIA thesis project</strong>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/R-analysis%20workflow-276DC3?logo=r&logoColor=white" alt="R analysis workflow">
  <img src="https://img.shields.io/badge/status-thesis%20analysis-6A5ACD" alt="Thesis analysis status">
  <img src="https://img.shields.io/badge/data-not%20included-important" alt="Data not included">
</p>

---

## Overview

This repository contains the R scripts used to process, quality-control, normalise, impute, model, and visualise maternal proteomics data from the ENDIA study.

The analyses investigate protein abundance patterns across pregnancy, with particular focus on:

- Maternal type 1 diabetes status
- Trimester-specific protein trajectories
- Maternal type 1 diabetes-by-trimester interactions
- Offspring case/control status
- Post-hoc estimated marginal means and multiple-testing correction

> **Data privacy note:** Individual-level proteomics data, participant metadata, and other protected study files are not included in this repository.

---

## Analysis workflow

```text
Raw proteomics data and protected metadata
                    ↓
        Data merging and protein filtering
                    ↓
 Missingness assessment and sample-level quality control
                    ↓
      Normalisation and RLE diagnostic assessment
                    ↓
           Missing-value imputation (QRILC)
                    ↓
  Mixed-effects modelling and multiple-testing correction
                    ↓
 Post-hoc contrasts, functional annotation, and visualisation
```

---

## Repository structure

```text
Proteomics-ENDIA/
│
├── R/
│   ├── 01_data_qc/           # Data merging, filtering, missingness and QC
│   ├── 02_preprocessing/     # Normalisation and imputation
│   ├── 03_models/            # Mixed-effects models, q-values and estimated means
│   ├── 04_visualization/     # Heatmaps and volcano plots
│   └── functions/            # Reusable multiple-testing helper functions
│
├── README.md                  # Project documentation
├── .gitignore                 # Excludes sensitive data and generated outputs
└── Proteomics-ENDIA.Rproj     # RStudio project file
```

### Main folders

| Folder | Purpose |
|---|---|
| `R/01_data_qc/` | Data merging, protein filtering, missing-value assessment, outlier identification, and QC figures |
| `R/02_preprocessing/` | Normalisation diagnostics, relative log expression plots, and QRILC imputation |
| `R/03_models/` | Mixed-effects models, false-discovery-rate correction, post-hoc contrasts, and estimated marginal means |
| `R/04_visualization/` | Heatmaps, volcano plots, and other results visualisations |
| `R/functions/` | Reusable q-value and multiple-testing functions |

---

## Reproducing the analysis

1. Clone this repository.

2. Open `Proteomics-ENDIA.Rproj` in RStudio.

3. Obtain approved access to the protected ENDIA data and metadata.

4. Configure the local file paths in the scripts. Protected data should remain outside this repository.

5. Run the workflow sequentially:

```text
R/01_data_qc/
        ↓
R/02_preprocessing/
        ↓
R/03_models/
        ↓
R/04_visualization/
```

---

## Statistical approach

The analysis workflow includes:

- Protein filtering based on presence thresholds
- Sample-level missingness and total-intensity quality control
- Outlier identification using multiple quality-control measures
- Variance-stabilising normalisation and RLE diagnostics
- QRILC imputation for missing proteomics values
- Linear mixed-effects modelling to account for repeated maternal samples
- False-discovery-rate correction for multiple testing
- Estimated marginal means and trimester-specific post-hoc contrasts

---

## Software

Analyses were conducted in R. Key packages include:

- `dplyr`
- `ggplot2`
- `lme4`
- `emmeans`
- `car`
- `vsn`
- `imputeLCMD`

---
