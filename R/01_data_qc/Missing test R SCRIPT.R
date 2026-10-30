library(dplyr)
library(tidyr)
library(ggplot2)
library(naniar)
library(stringr)
library(tibble)


#-------------------------------------------------------------------
# Assumes `unique_samples`, `intensity_mat`, `missing_summary`, and
# `prot_cols` already exist from the previous QC script (now sourced
# from protein_column_names.rds, not regex — see that script's fix)
#-------------------------------------------------------------------

# Rebuild missing_summary on the FILTERED (451-protein) matrix —
# the old missing_summary was computed on all 469 proteins, pre-filter
missing_summary <- data.frame(
  Protein     = colnames(intensity_mat),
  n_missing   = colSums(is.na(intensity_mat)),
  pct_missing = colMeans(is.na(intensity_mat)) * 100
)

#intensity_mat <- intensity_mat_filtered   # so the rest of the script needs no other edits

#===================================================================
# 1. Little's MCAR test (overall omnibus check)
#===================================================================

# naniar::mcar_test needs a data frame of numeric variables (samples x proteins)
#mcar_result <- mcar_test(as.data.frame(intensity_mat))
#print(mcar_result)
# p < 0.05 => reject MCAR (missingness is NOT completely random overall)
# Note: with many variables/proteins this test is very sensitive and
# will often reject MCAR even with fairly mild, patterned missingness —
# treat it as one piece of evidence, not the final word.

#===================================================================
# 2. MNAR signature: missingness vs. protein abundance
#===================================================================

# Mean intensity of each protein, computed only from its OBSERVED
# (non-missing) values
mean_intensity <- colMeans(intensity_mat, na.rm = TRUE)

mnar_check <- missing_summary %>%
  mutate(mean_log2_intensity = log2(mean_intensity[Protein] + 1))

# Correlation: negative => lower-abundance proteins go missing more often (MNAR signature)
cor_test <- cor.test(mnar_check$pct_missing, mnar_check$mean_log2_intensity,
                     method = "spearman")
print(cor_test)

ggplot(mnar_check, aes(x = mean_log2_intensity, y = pct_missing)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, colour = "firebrick") +
  theme_minimal() +
  labs(title = "Missingness vs. mean protein abundance",
       subtitle = paste0("Spearman rho = ", round(cor_test$estimate, 2),
                         ", p = ", signif(cor_test$p.value, 3)),
       x = "Mean log2 intensity (observed values only)",
       y = "% missing")

ggsave("missingness_vs_abundance.png", width = 7, height = 5)

#===================================================================
# 3. MAR signature: missingness vs. known covariates
#===================================================================

# *** FIX ***
# Plate.x was being fed into glm() as-is. If it's stored as an
# integer/numeric column (1, 2, 3), glm() silently fits it as a
# CONTINUOUS predictor — i.e. one linear slope across plate number —
# rather than estimating a separate effect per plate. That's not what
# your writeup reports ("Plate 2 — significant, p = 1.55e-7"), which
# implies per-plate levels were actually estimated. Cast it to a
# factor explicitly so each plate gets its own coefficient.
meta_small <- unique_samples %>%
  select(SampleID, case, Plate.x, corrected_visit_sample) %>%
  mutate(Plate.x = factor(Plate.x))

miss_long <- as.data.frame(is.na(intensity_mat)) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Protein", values_to = "is_missing") %>%
  left_join(meta_small, by = "SampleID")

# Overall missingness rate by group — eyeball first
miss_long %>%
  group_by(case) %>%
  summarise(pct_missing = mean(is_missing) * 100)

miss_long %>%
  group_by(Plate.x) %>%
  summarise(pct_missing = mean(is_missing) * 100)

miss_long %>%
  group_by(corrected_visit_sample) %>%
  summarise(pct_missing = mean(is_missing) * 100)

# Formal check: logistic regression of missingness ~ covariates,
# fit ACROSS all proteins together (quick omnibus check).
# Plate.x is now a factor (see fix above), so this gives one
# coefficient per plate rather than a single linear-trend coefficient.
mar_model <- glm(is_missing ~ case + Plate.x + corrected_visit_sample,
                 data = miss_long, family = binomial)
summary(mar_model)
# Significant coefficients => missingness is associated with that
# covariate independent of protein identity => evidence for MAR
# driven by that variable (e.g. Plate.x/batch effects)

#===================================================================
# Overall missingness — single summary number
#===================================================================

overall_pct_missing <- mean(is.na(intensity_mat)) * 100

cat("Overall missingness:", round(overall_pct_missing, 2), "%\n")
cat("Missing values:", sum(is.na(intensity_mat)), "of",
    length(intensity_mat), "total data points\n")

#===================================================================
# 4. Per-protein view: which proteins look most MNAR-like?
#===================================================================

# Rank proteins by "MNAR-ness": high missingness + low observed abundance
mnar_check <- mnar_check %>%
  arrange(desc(pct_missing), mean_log2_intensity)

cat("\nProteins most consistent with MNAR (high missing, low abundance):\n")
print(head(mnar_check %>% select(Protein, pct_missing, mean_log2_intensity), 15))