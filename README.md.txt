# ENDIA Proteomics Analysis

> Reproducible R workflows for the maternal proteomics component of the ENDIA thesis project.

![R](https://img.shields.io/badge/R-analysis%20workflow-276DC3)
![Status](https://img.shields.io/badge/status-thesis%20analysis-6A5ACD)
![Data](https://img.shields.io/badge/data-not%20included-important)

## Overview

This repository contains the R analysis workflows used to process, quality-control, normalise, impute, model, and visualise maternal proteomics data from the ENDIA study.

The project investigates protein abundance patterns across pregnancy, including maternal type 1 diabetes status, trimester effects, and offspring case/control outcomes.

> Individual-level data and identifiable metadata are not included in this repository.

## Analysis workflow

```text
Raw proteomics and metadata
          ↓
Data merging and protein filtering
          ↓
Missingness assessment and sample quality control
          ↓
Normalisation and relative log expression diagnostics
          ↓
Missing-value imputation
          ↓
Mixed-effects modelling and multiple-testing correction
          ↓
Post-hoc contrasts, functional annotation, and visualisation