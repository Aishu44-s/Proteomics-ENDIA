# ============================================================================
# D6_bulk_emmeans_t1d_visit.R
# Adjusted trimester trajectories and pairwise trimester contrasts
# for proteins significant for the T1D-by-trimester interaction
# ============================================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

# ---------------------------------------------------------------------------
# Choose which significant result set to plot
# ---------------------------------------------------------------------------

target_term <- "corrected_visit_sample:maternal_T1D"
q_threshold <- 0.05

# Read the master results table created by D5
fixed_effects_master <- read.csv(
  "fixed_effects_master_results_adj_int_rm.csv",
  check.names = FALSE
)

# Select significant proteins for the chosen interaction
sig_proteins <- fixed_effects_master %>%
  filter(
    Variable == target_term,
    QValue < q_threshold
  ) %>%
  pull(Protein) %>%
  unique()

length(sig_proteins)
sig_proteins

if (length(sig_proteins) == 0) {
  stop("No proteins meet the selected q-value threshold.")
}

# ---------------------------------------------------------------------------
# Use exactly the same adjusted model as your D3/D5 result table
# ---------------------------------------------------------------------------

model_name <- "form_adj_int_rm"
model_formula <- as.formula(model[[model_name]])

# ---------------------------------------------------------------------------
# Empty lists to collect results
# ---------------------------------------------------------------------------

emm_means_list <- list()
emm_contrasts_list <- list()
fit_log <- list()

# ---------------------------------------------------------------------------
# Fit one model per significant protein and extract emmeans
# ---------------------------------------------------------------------------

for (protein in sig_proteins) {
  
  cat("Processing:", protein, "\n")
  
  # Make a temporary dataset with the selected protein as btmean
  tmp <- mother
  tmp$btmean <- tmp[[protein]]
  tmp <- tmp[!is.na(tmp$btmean), ]
  
  # Fit the same adjusted model used for the master table
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
  
  # Record and skip failed fits
  if (inherits(fit, "error")) {
    
    fit_log[[protein]] <- data.frame(
      Protein = protein,
      Status = "Failed",
      Message = fit$message,
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  fit_log[[protein]] <- data.frame(
    Protein = protein,
    Status = "Success",
    Message = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # -------------------------------------------------------------------------
  # Adjusted mean protein abundance:
  # One mean per trimester, separately for each maternal-T1D group
  # -------------------------------------------------------------------------
  
  emm <- emmeans(
    fit,
    ~ corrected_visit_sample | maternal_T1D,
    lmer.df = "satterthwaite"
  )
  
  means_df <- as.data.frame(
    summary(emm, infer = c(TRUE, TRUE))
  )
  
  means_df$Protein <- protein
  means_df$ModelName <- model_name
  
  emm_means_list[[protein]] <- means_df
  
