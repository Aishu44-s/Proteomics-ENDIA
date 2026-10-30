# D5: Create fixed-effect results tables after running D3

library(dplyr)
library(qvalue)

if (!file.exists("qvaluetruncfunction.R")) {
  stop("Cannot find qvaluetruncfunction.R in the working directory.")
}
source("qvaluetruncfunction.R")

# Confirm D3 has been run in this R session
required_objects <- c("model", "resp", "anova.ls")
missing_objects <- required_objects[
  !vapply(required_objects, exists, logical(1), inherits = TRUE)
]

if (length(missing_objects) > 0) {
  stop(
    "Run D3_analysis_mother_btmean.R first in this R session. Missing: ",
    paste(missing_objects, collapse = ", ")
  )
}

# Choose the adjusted model with the three-way interaction
model_name <- "form_adj_int_rm"
model_index <- match(model_name, names(model))

if (is.na(model_index)) {
  stop("Could not find ", model_name, " in the model list.")
}

# Primary term for the research question
primary_term <- "casecat:corrected_visit_sample:maternal_T1D"

# Terms to test. The three-way interaction is primary;
# the other terms are secondary.
terms_to_test <- c(
  "casecat",
  "corrected_visit_sample",
  "maternal_T1D",
  "casecat:corrected_visit_sample",
  "casecat:maternal_T1D",
  "corrected_visit_sample:maternal_T1D",
  primary_term,
  "infant_sex",
  "timefromonsetmonths",
  "maternal_HLA",
  "maternal_age_at_conception",
  "maternal_conception_bmi",
  "maternal_smoking"
)

# Find an ANOVA table from a successfully fitted protein
available_terms <- NULL

for (i in seq_along(resp)) {
  anova_object <- anova.ls[[i]][[model_index]]
  
  if (!is.null(anova_object)) {
    available_terms <- rownames(as.data.frame(anova_object))
    break
  }
}

if (is.null(available_terms)) {
  stop("No ANOVA results were found for model ", model_name, ".")
}

missing_terms <- setdiff(terms_to_test, available_terms)

if (length(missing_terms) > 0) {
  stop(
    "These terms were not found in the ANOVA results: ",
    paste(missing_terms, collapse = ", "),
    "\nAvailable terms: ",
    paste(available_terms, collapse = ", ")
  )
}

# Extract one p-value for a term from an ANOVA table
get_anova_p <- function(anova_object, term) {
  if (is.null(anova_object)) {
    return(NA_real_)
  }
  
  anova_table <- as.data.frame(anova_object)
  p_column <- grep("^Pr\\(", names(anova_table), value = TRUE)
  
  if (length(p_column) != 1) {
    return(NA_real_)
  }
  
  row_index <- which(rownames(anova_table) == term)
  
  if (length(row_index) != 1) {
    return(NA_real_)
  }
  
  as.numeric(anova_table[row_index, p_column])
}

# Calculate q-values for one term across all proteins
safe_qvalues <- function(p) {
  q_out <- rep(NA_real_, length(p))
  valid_p <- is.finite(p) & p >= 0 & p <= 1
  
  if (sum(valid_p) == 0) {
    return(list(
      qvalues = q_out,
      method = "No valid p-values"
    ))
  }
  
  q_attempts <- list(
    qvalue_default = function(x) {
      qvalue_truncp(x, lfdr.out = FALSE)
    },
    qvalue_lambda_0.05_to_0.50 = function(x) {
      qvalue_truncp(
        x,
        lambda = seq(0.05, 0.50, 0.05),
        lfdr.out = FALSE
      )
    },
    qvalue_lambda_0.50 = function(x) {
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
  
  # Use Benjamini-Hochberg if the q-value methods fail
  q_out[valid_p] <- p.adjust(p[valid_p], method = "BH")
  
  list(
    qvalues = q_out,
    method = "BH_fallback"
  )
}

# Collect raw p-values: rows are terms and columns are proteins
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

# Calculate q-values separately for each term across proteins
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

# Create one row per protein and tested term
fixed_effects_master <- bind_rows(
  lapply(terms_to_test, function(term) {
    data.frame(
      Variable = term,
      ResultType = ifelse(term == primary_term, "Primary", "Secondary"),
      Protein = resp,
      ModelName = model_name,
      Formula = unname(model[[model_index]]),
      RawPValue = as.numeric(pvalue_mat[term, ]),
      QValue = as.numeric(qvalue_mat[term, ]),
      Significant_q_0_10 = as.numeric(qvalue_mat[term, ]) < 0.10,
      QValueMethod = unname(qvalue_method[term]),
      stringsAsFactors = FALSE
    )
  })
) %>%
  arrange(ResultType, Variable, QValue)

# Summarise the number of proteins tested and significant for each term
fixed_effects_hit_table <- fixed_effects_master %>%
  group_by(ResultType, Variable, ModelName, QValueMethod) %>%
  summarise(
    ProteinsTested = sum(!is.na(RawPValue)),
    Q_lt_0_10 = sum(Significant_q_0_10, na.rm = TRUE),
    .groups = "drop"
  )

# Primary three-way interaction results, including all proteins
primary_interaction_results <- fixed_effects_master %>%
  filter(Variable == primary_term) %>%
  arrange(QValue)

primary_interaction_hits_q10 <- primary_interaction_results %>%
  filter(Significant_q_0_10)

# Secondary results for the other terms
secondary_results <- fixed_effects_master %>%
  filter(Variable != primary_term)

# Save results
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
  primary_interaction_results,
  "primary_three_way_interaction_results_adj_int_rm.csv",
  row.names = FALSE
)

write.csv(
  primary_interaction_hits_q10,
  "primary_three_way_interaction_significant_q10_adj_int_rm.csv",
  row.names = FALSE
)

write.csv(
  secondary_results,
  "secondary_fixed_effects_results_adj_int_rm.csv",
  row.names = FALSE
)

message("D5 results tables saved successfully.")
message("Primary term: ", primary_term)