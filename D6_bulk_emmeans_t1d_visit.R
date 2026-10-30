# ============================================================================
# D6_bulk_emmeans_t1d_visit.R
# Run after D3/D5, in an R session containing mother and model
# ============================================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

target_term <- "corrected_visit_sample:maternal_T1D"
q_threshold <- 0.10

model_name <- "form_adj_int_rm"
pairwise_adjustment <- "tukey"

# ---------------------------------------------------------------------------
# Read D5 master results and select significant proteins
# ---------------------------------------------------------------------------

fixed_effects_master <- read.csv(
  "fixed_effects_master_results_adj_int_rm.csv",
  check.names = FALSE
)

selected_proteins <- fixed_effects_master %>%
  filter(
    Variable == target_term,
    QValue < q_threshold
  ) %>%
  select(Protein, QValue) %>%
  distinct() %>%
  arrange(QValue)

sig_proteins <- selected_proteins$Protein

if (length(sig_proteins) == 0) {
  stop("No proteins meet the selected q-value threshold.")
}

message("Number of selected proteins: ", length(sig_proteins))

# ---------------------------------------------------------------------------
# Use the same adjusted model as D3/D5
# ---------------------------------------------------------------------------

if (!model_name %in% names(model)) {
  stop(paste("Could not find", model_name, "in model."))
}

model_formula <- as.formula(model[[model_name]])

# ---------------------------------------------------------------------------
# Empty lists for results
# ---------------------------------------------------------------------------

emm_means_list <- list()
emm_contrasts_list <- list()
fit_log_list <- list()
design_counts_list <- list()
nonestimable_list <- list()

# ---------------------------------------------------------------------------
# Fit one model per protein and extract estimated marginal means
# ---------------------------------------------------------------------------

for (protein in sig_proteins) {
  
  cat("Processing:", protein, "\n")
  
  protein_q <- selected_proteins %>%
    filter(Protein == protein) %>%
    pull(QValue)
  
  # Create temporary data with the protein as the outcome
  tmp <- mother
  tmp$btmean <- tmp[[protein]]
  tmp <- tmp[!is.na(tmp$btmean), ]
  
  # Fit the same adjusted model used in D3
  fit <- tryCatch(
    lmer(
      formula = model_formula,
      data = tmp,
      weights = weightTEDDY,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) e
  )
  
  # Save model errors and move to the next protein
  if (inherits(fit, "error")) {
    
    fit_log_list[[protein]] <- data.frame(
      Protein = protein,
      InteractionQValue = protein_q,
      Status = "Failed",
      Message = fit$message,
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  fit_log_list[[protein]] <- data.frame(
    Protein = protein,
    InteractionQValue = protein_q,
    Status = "Success",
    Message = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Diagnostic: counts actually used by this fitted model
  # This uses model.frame(fit), not the original mother dataset.
  # -------------------------------------------------------------------------
  
  model_data <- model.frame(fit)
  
  design_counts <- as.data.frame(
    with(
      model_data,
      table(
        casecat,
        maternal_T1D,
        corrected_visit_sample,
        useNA = "ifany"
      )
    )
  )
  
  names(design_counts)[4] <- "N"
  
  design_counts$Protein <- protein
  design_counts$InteractionQValue <- protein_q
  
  design_counts_list[[protein]] <- design_counts
  
  # -------------------------------------------------------------------------
  # Estimated means by trimester within maternal-T1D group
  # -------------------------------------------------------------------------
  
  analysis_casecat <- "0"  # Confirm what 0 means before presenting
  
  emm <- tryCatch(
    emmeans(
      fit,
      ~ corrected_visit_sample | maternal_T1D,
      at = list(casecat = analysis_casecat),
      lmer.df = "satterthwaite"
    ),
    error = function(e) e
  )
  
  if (inherits(emm, "error")) {
    
    fit_log_list[[protein]] <- data.frame(
      Protein = protein,
      InteractionQValue = protein_q,
      Status = "emmeans_failed",
      Message = emm$message,
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  means_df <- as.data.frame(
    summary(
      emm,
      infer = c(TRUE, TRUE),
      level = 0.90
    )
  )
  
  means_df$Protein <- protein
  means_df$ModelName <- model_name
  means_df$InteractionQValue <- protein_q
  means_df$casecat_condition <- analysis_casecat
  
  emm_means_list[[protein]] <- means_df
  
  # Record cells where emmeans could not calculate an adjusted mean
  nonestimable_rows <- means_df %>%
    filter(is.na(emmean)) %>%
    select(
      Protein,
      ModelName,
      InteractionQValue,
      maternal_T1D,
      corrected_visit_sample
    )
  
  if (nrow(nonestimable_rows) > 0) {
    nonestimable_list[[protein]] <- nonestimable_rows
  }
  
  # -------------------------------------------------------------------------
  # All trimester-pair comparisons within each maternal-T1D group
  #
  # With 3 trimesters this includes:
  # T1 vs T2, T1 vs T3, and T2 vs T3
  # separately for each maternal-T1D group.
  # -------------------------------------------------------------------------
  
  contrasts_df <- as.data.frame(
    summary(
      pairs(emm, adjust = pairwise_adjustment),
      infer = c(TRUE, TRUE),
      level = 0.90
    )
  )
  
  contrasts_df$Protein <- protein
  contrasts_df$ModelName <- model_name
  contrasts_df$InteractionQValue <- protein_q
  contrasts_df$casecat_condition <- analysis_casecat
  
  emm_contrasts_list[[protein]] <- contrasts_df
}

# ---------------------------------------------------------------------------
# Combine all results
# ---------------------------------------------------------------------------

emm_means_long <- bind_rows(emm_means_list)

emm_contrasts_long <- bind_rows(emm_contrasts_list)

emm_fit_log <- bind_rows(fit_log_list)

emm_design_counts <- bind_rows(design_counts_list)

nonestimable_emmeans <- bind_rows(nonestimable_list)

# ---------------------------------------------------------------------------
# Save output files
# ---------------------------------------------------------------------------

write.csv(
  emm_means_long,
  "emmeans_t1d_visit_casecat0_means_q10.csv",
  row.names = FALSE
)

write.csv(
  emm_contrasts_long,
  "emmeans_t1d_visit_casecat0_pairwise_contrasts_q10.csv",
  row.names = FALSE
)

write.csv(
  emm_fit_log,
  "emmeans_t1d_visit_fit_log_q10.csv",
  row.names = FALSE
)

write.csv(
  emm_design_counts,
  "emmeans_t1d_visit_design_counts_q10.csv",
  row.names = FALSE
)

emm_design_counts_casecat0 <- emm_design_counts %>%
  filter(casecat == analysis_casecat)

write.csv(
  emm_design_counts_casecat0,
  "emmeans_t1d_visit_casecat0_design_counts_q10.csv",
  row.names = FALSE
)

write.csv(
  nonestimable_emmeans,
  "emmeans_t1d_visit_casecat0_nonestimable_cells_q10.csv",
  row.names = FALSE
)

message("Bulk emmeans extraction for casecat = 0 completed.")
message("Check emmeans_t1d_visit_casecat0_design_counts_q10.csv for counts.")
message("Check emmeans_t1d_visit_casecat0_nonestimable_cells_q10.csv for unavailable adjusted means.")