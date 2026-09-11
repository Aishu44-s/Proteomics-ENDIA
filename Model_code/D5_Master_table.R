# ============================================================================
# D5_fixed_effects_significance.R
# Run AFTER D3_analysis_mother_btmean.R has completed
#
# Creates a filterable Excel-friendly table for:
# - case/control
# - trimester / visit
# - maternal T1D
# - all biological interactions
# - confounders
#
# Multiple-testing correction:
# q-values are calculated separately for each term across all proteins.
# ============================================================================

library(dplyr)
library(qvalue)

source("qvaluetruncfunction.R")

# ---------------------------------------------------------------------------
# Choose the model to report
# ---------------------------------------------------------------------------

model_name <- "form_adj_int_rm"

model_index <- match(model_name, names(model))

if (is.na(model_index)) {
  stop(paste("Could not find", model_name, "in names(model)."))
}

# ---------------------------------------------------------------------------
# Terms to test
#
# These must exactly match the row names in car::Anova() output.
# ---------------------------------------------------------------------------

terms_to_test <- c(
  
  # Main biological effects
  "casecat",
  "corrected_visit_sample",
  "maternal_T1D",
  
  # Biological interactions
  "casecat:corrected_visit_sample",
  "casecat:maternal_T1D",
  "corrected_visit_sample:maternal_T1D",
  "casecat:corrected_visit_sample:maternal_T1D",
  
  # Confounders
  "infant_sex",
  "timefromonsetmonths",
  "maternal_HLA",
  "maternal_age_at_conception",
  "maternal_conception_bmi",
  "maternal_smoking"
)

# Check the exact available ANOVA term labels before proceeding
available_terms <- rownames(as.data.frame(anova.ls[[1]][[model_index]]))

missing_terms <- setdiff(terms_to_test, available_terms)

if (length(missing_terms) > 0) {
  stop(
    paste(
      "These terms were not found in the ANOVA table:",
      paste(missing_terms, collapse = ", "),
      "\nAvailable terms are:",
      paste(available_terms, collapse = ", ")
    )
  )
}

# ---------------------------------------------------------------------------
# Helper: extract one Type II Wald ANOVA p-value
# ---------------------------------------------------------------------------

get_anova_p <- function(anova_object, term) {
  
  if (is.null(anova_object)) return(NA_real_)
  
  anova_table <- as.data.frame(anova_object)
  p_column <- grep("^Pr\\(", names(anova_table), value = TRUE)
  
  if (length(p_column) != 1) return(NA_real_)
  
  row_index <- which(rownames(anova_table) == term)
  
  if (length(row_index) != 1) return(NA_real_)
  
  as.numeric(anova_table[row_index, p_column])
}

# ---------------------------------------------------------------------------
# Helper: calculate q-values robustly
# ---------------------------------------------------------------------------

safe_qvalues <- function(p) {
  
  q_out <- rep(NA_real_, length(p))
  valid_p <- !is.na(p)
  
  if (sum(valid_p) == 0) {
    return(list(
      qvalues = q_out,
      method = "No valid p-values"
    ))
  }
  
  q_attempts <- list(
    
    "qvalue_default" = function(x) {
      qvalue_truncp(x, lfdr.out = FALSE)
    },
    
    "qvalue_lambda_0.05_to_0.50" = function(x) {
      qvalue_truncp(
        x,
        lambda = seq(0.05, 0.50, 0.05),
        lfdr.out = FALSE
      )
    },
    
    "qvalue_lambda_0.50" = function(x) {
      qvalue_truncp(
        x,
        lambda = 0.50,
        lfdr.out = FALSE
      )
    }
  )
  
  for (method_name in names(q_attempts)) {
    
    q_result <- tryCatch(
      q_attempts[[method_name]](p[valid_p]),
      error = function(e) NULL
    )
    
    if (!is.null(q_result)) {
      
      q_out[valid_p] <- q_result$qvalues
      
      return(list(
        qvalues = q_out,
        method = method_name
      ))
    }
  }
  
  # Only used if every qvalue method fails.
  q_out[valid_p] <- p.adjust(p[valid_p], method = "BH")
  
  list(
    qvalues = q_out,
    method = "BH_fallback"
  )
}

# ---------------------------------------------------------------------------
# Extract raw p-values
# Rows = model terms; columns = proteins
# ---------------------------------------------------------------------------

pvalue_mat <- matrix(
  NA_real_,
  nrow = length(terms_to_test),
  ncol = length(resp),
  dimnames = list(terms_to_test, resp)
)

for (term in terms_to_test) {
  
  for (i in seq_along(resp)) {
    
    pvalue_mat[term, i] <- get_anova_p(
      anova_object = anova.ls[[i]][[model_index]],
      term = term
    )
  }
}

# ---------------------------------------------------------------------------
# Calculate q-values separately for each term across all proteins
# ---------------------------------------------------------------------------

qvalue_mat <- matrix(
  NA_real_,
  nrow = length(terms_to_test),
  ncol = length(resp),
  dimnames = list(terms_to_test, resp)
)

qvalue_method <- setNames(
  rep(NA_character_, length(terms_to_test)),
  terms_to_test
)

for (term in terms_to_test) {
  
  q_result <- safe_qvalues(pvalue_mat[term, ])
  
  qvalue_mat[term, ] <- q_result$qvalues
  qvalue_method[term] <- q_result$method
}

# ---------------------------------------------------------------------------
# Create the main Excel-friendly master table
# One row = one term x one protein
# ---------------------------------------------------------------------------

fixed_effects_master <- bind_rows(
  
  lapply(terms_to_test, function(term) {
    
    data.frame(
      Variable = term,
      Protein = resp,
      ModelName = model_name,
      Formula = unname(model[[model_index]]),
      RawPValue = as.numeric(pvalue_mat[term, ]),
      QValue = as.numeric(qvalue_mat[term, ]),
      RawP_lt_0_05 = as.numeric(pvalue_mat[term, ]) < 0.05,
      Significant_q_0_05 = as.numeric(qvalue_mat[term, ]) < 0.05,
      Significant_q_0_10 = as.numeric(qvalue_mat[term, ]) < 0.10,
      QValueMethod = qvalue_method[term],
      stringsAsFactors = FALSE
    )
  })
)

# Sort so significant proteins appear first for each variable
fixed_effects_master <- fixed_effects_master %>%
  arrange(Variable, QValue)

# ---------------------------------------------------------------------------
# Hit-count summary: one row per variable
# ---------------------------------------------------------------------------

fixed_effects_hit_table <- fixed_effects_master %>%
  group_by(Variable, ModelName, QValueMethod) %>%
  summarise(
    ProteinsTested = sum(!is.na(RawPValue)),
    RawP_lt_0_05 = sum(RawP_lt_0_05, na.rm = TRUE),
    Q_lt_0_05 = sum(Significant_q_0_05, na.rm = TRUE),
    Q_lt_0_10 = sum(Significant_q_0_10, na.rm = TRUE),
    .groups = "drop"
  )

# Already-filtered significant-protein tables
fixed_effects_hits_q05 <- fixed_effects_master %>%
  filter(Significant_q_0_05)

fixed_effects_hits_q10 <- fixed_effects_master %>%
  filter(Significant_q_0_10)

# ---------------------------------------------------------------------------
# Save files
# ---------------------------------------------------------------------------

write.csv(
  pvalue_mat,
  "fixed_effects_pvalues_adj_int_rm.csv"
)

write.csv(
  qvalue_mat,
  "fixed_effects_qvalues_adj_int_rm.csv"
)

write.csv(
  fixed_effects_master,
  "fixed_effects_master_results_adj_int_rm.csv",
  row.names = FALSE
)

write.csv(
  fixed_effects_hit_table,
  "fixed_effects_hit_table_adj_int_rm.csv",
  row.names = FALSE
)

write.csv(
  fixed_effects_hits_q05,
  "fixed_effects_significant_q05_adj_int_rm.csv",
  row.names = FALSE
)

write.csv(
  fixed_effects_hits_q10,
  "fixed_effects_significant_q10_adj_int_rm.csv",
  row.names = FALSE
)

message("Fixed-effect master tables saved successfully.")