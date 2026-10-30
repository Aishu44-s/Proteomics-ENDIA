
library(tidyverse)
library(pcaMethods)  
library(ggrepel)
library(patchwork)

#-------------------------------------------------------------------
# 1. Load data — pre-imputation, still has NAs
#-------------------------------------------------------------------

protein_filtered <- read.csv("protein_filtered_80pct.csv", header = TRUE)
prot_cols_filt   <- readRDS("protein_column_names_post_filter.rds")

 meta_clean <- protein_filtered %>%
  select(-all_of(prot_cols_filt)) %>%
  rename(Plate = Plate.x)

protein_log2_clean <- protein_filtered %>%
  select(all_of(prot_cols_filt)) %>%
  as.matrix() %>%
  {log2(. + 1)}
rownames(protein_log2_clean) <- protein_filtered$SampleID

cat("Samples (pre-imputation):", nrow(protein_log2_clean),
    "| Proteins:", ncol(protein_log2_clean),
    "| NAs:", sum(is.na(protein_log2_clean)), "\n")

#-------------------------------------------------------------------
# 2. NIPALS-PCA (handles NAs natively — no imputation needed first)
#-------------------------------------------------------------------

pca_res <- pca(protein_log2_clean, method = "nipals", nPcs = 2,
               scale = "uv", center = TRUE)

pca_df <- as.data.frame(scores(pca_res)) %>%
  rownames_to_column("SampleID") %>%
  left_join(
    meta_clean %>%
      select(SampleID, mother_id, Plate, case, corrected_visit_sample,
             maternal_T1D, maternal_age_at_conception,
             maternal_conception_bmi, ga_at_birth_weeks),
    by = "SampleID"
  ) %>%
  mutate(
    pc1_z = as.numeric(scale(PC1)),
    pc2_z = as.numeric(scale(PC2)),
    outlier_zscore = abs(pc1_z) > 3 | abs(pc2_z) > 3
  )

var1 <- round(pca_res@R2[1] * 100, 1)
var2 <- round(pca_res@R2[2] * 100, 1)

#-------------------------------------------------------------------
# 3. Multi-method outlier flagging
#-------------------------------------------------------------------

# -- Mahalanobis distance (multivariate, on PC1+PC2) --
pc_mat <- as.matrix(pca_df[, c("PC1", "PC2")])
pca_df$mahal_dist   <- mahalanobis(pc_mat, colMeans(pc_mat), cov(pc_mat))
pca_df$mahal_pval   <- pchisq(pca_df$mahal_dist, df = 2, lower.tail = FALSE)
pca_df$outlier_mahal <- pca_df$mahal_pval < 0.001

# -- Total intensity, computed on PRE-imputation data (na.rm = TRUE,
#    i.e. sum of observed values only — consistent with using
#    pre-imputation data throughout this script) --
pca_df <- pca_df %>%
  mutate(
    total_intensity  = rowSums(protein_log2_clean[SampleID, ], na.rm = TRUE),
    intensity_z       = as.numeric(scale(total_intensity)),
    outlier_intensity = intensity_z < -3
  )

# -- % missing (>50% threshold) --
miss_rate <- rowMeans(is.na(protein_log2_clean)) * 100
pca_df$pct_missing     <- miss_rate[pca_df$SampleID]
pca_df$outlier_missing <- pca_df$pct_missing > 50

# -- IQR method on PC1 --
Q1 <- quantile(pca_df$PC1, 0.25); Q3 <- quantile(pca_df$PC1, 0.75)
IQR_val <- Q3 - Q1
pca_df$outlier_iqr <- pca_df$PC1 < (Q1 - 3 * IQR_val) | pca_df$PC1 > (Q3 + 3 * IQR_val)

pca_df <- pca_df %>%
  mutate(
    n_methods_flagged = outlier_zscore + outlier_mahal + outlier_intensity +
      outlier_missing + outlier_iqr
  )

strong_outliers <- pca_df %>%
  filter(n_methods_flagged >= 2) %>%
  arrange(desc(n_methods_flagged), desc(mahal_dist)) %>%
  select(SampleID, mother_id, Plate, case, corrected_visit_sample,
         pc1_z, pc2_z, mahal_dist, mahal_pval,
         pct_missing, total_intensity, n_methods_flagged)

cat("\nSamples flagged by >=2 methods:\n")
print(strong_outliers)
write.csv(strong_outliers, "outlier_investigation_summary.csv", row.names = FALSE)

#-------------------------------------------------------------------
# 4. Diagnostic plots
#-------------------------------------------------------------------

p_pca_labelled <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = n_methods_flagged >= 2,
                                     shape = factor(case))) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_text_repel(data = filter(pca_df, n_methods_flagged >= 2),
                  aes(label = SampleID), size = 3, max.overlaps = 20, colour = "black") +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
                      labels = c("Normal", "Outlier (>=2 methods)")) +
  scale_shape_manual(values = c("0" = 16, "1" = 17), labels = c("Control", "Case")) +
  labs(title = "NIPALS-PCA (pre-imputation) — Outlier Samples Labelled",
       colour = "Status", shape = "Group",
       x = paste0("PC1 (", var1, "%)"), y = paste0("PC2 (", var2, "%)")) +
  theme_classic(base_size = 12)
ggsave("pca_outliers_labelled.png", p_pca_labelled, width = 9, height = 7, dpi = 600)

p1 <- ggplot(pca_df, aes(x = reorder(SampleID, mahal_dist), y = mahal_dist, colour = outlier_mahal)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = qchisq(0.999, df = 2), linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  labs(title = "Mahalanobis Distance (PC1+PC2)", subtitle = "Red line = p<0.001 threshold",
       x = NULL, y = "Mahalanobis Distance") +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none")

p2 <- ggplot(pca_df, aes(x = total_intensity, fill = outlier_intensity)) +
  geom_histogram(bins = 40, alpha = 0.8, colour = "white") +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  labs(title = "Total Intensity (observed values only)", subtitle = "Red = >3 SD below mean",
       x = "Total log2 Intensity", y = "Count") +
  theme_classic(base_size = 11) + theme(legend.position = "none")

p3 <- ggplot(pca_df, aes(x = pct_missing, fill = outlier_missing)) +
  geom_histogram(bins = 30, alpha = 0.8, colour = "white") +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick")) +
  labs(title = "% Missing Values per Sample", subtitle = "Red line = 50% threshold",
       x = "% Missing", y = "Count") +
  theme_classic(base_size = 11) + theme(legend.position = "none")

p4 <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = factor(n_methods_flagged),
                         size = factor(n_methods_flagged))) +
  geom_point(alpha = 0.8) +
  geom_text_repel(data = filter(pca_df, n_methods_flagged >= 2),
                  aes(label = SampleID), size = 2.8, colour = "black", max.overlaps = 20) +
  scale_colour_manual(values = c("0" = "grey70", "1" = "orange", "2" = "red",
                                 "3" = "darkred", "4" = "purple", "5" = "black"),
                      name = "# Methods\nFlagged") +
  scale_size_manual(values = c("0"=1.5,"1"=2.5,"2"=3.5,"3"=4,"4"=4.5,"5"=5), guide = "none") +
  labs(title = "PCA — Outlier Convergence Across Methods",
       x = paste0("PC1 (", var1, "%)"), y = paste0("PC2 (", var2, "%)")) +
  theme_classic(base_size = 11)

combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Outlier Investigation — Multi-Method Analysis (NIPALS-PCA, pre-imputation)")
ggsave("outlier_multimethod_panel.png", combined, width = 14, height = 10, dpi = 650)

cat("Saved: pca_outliers_labelled.png, outlier_multimethod_panel.png\n")

#-------------------------------------------------------------------
# 5. Missing-protein pattern among flagged outliers
#-------------------------------------------------------------------

outlier_ids <- strong_outliers$SampleID
outlier_miss <- protein_log2_clean[rownames(protein_log2_clean) %in% outlier_ids, , drop = FALSE]

miss_long <- as.data.frame(is.na(outlier_miss)) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Protein", values_to = "Missing") %>%
  filter(Missing == TRUE) %>%
  mutate(short_name = gsub(".*\\.", "", Protein))

miss_summary <- miss_long %>%
  count(SampleID, name = "n_missing") %>%
  mutate(pct_missing = round(n_missing / ncol(outlier_miss) * 100, 1)) %>%
  arrange(desc(n_missing))
cat("\nMissing proteins per flagged outlier sample:\n")
print(miss_summary)

write.csv(miss_long %>% select(SampleID, short_name, Protein),
          "outlier_missing_proteins.csv", row.names = FALSE)

#-------------------------------------------------------------------
# 6. FINAL REMOVAL — all samples flagged by >=2 methods, per
#    supervisor: these are failed samples and should all be removed.
#-------------------------------------------------------------------

final_remove_ids <- strong_outliers %>%
  arrange(desc(n_methods_flagged), desc(mahal_dist)) %>%
  pull(SampleID)

cat("\n*** Samples marked for removal (all flagged by >=2 methods):", length(final_remove_ids), "***\n")
print(final_remove_ids)

#-------------------------------------------------------------------
# 7. Save the reduced dataset, ready for the imputation script to
#    consume directly — same format as protein_filtered_80pct.csv
#-------------------------------------------------------------------

protein_filtered_reduced <- protein_filtered %>%
  filter(!SampleID %in% final_remove_ids)

cat("\nSamples before removal:", nrow(protein_filtered), "\n")
cat("Samples after removal: ", nrow(protein_filtered_reduced), "\n")

write.csv(protein_filtered_reduced, "protein_filtered_80pct_outlier_removed.csv",
          row.names = FALSE)
saveRDS(final_remove_ids, "outlier_removed_sample_ids.rds")

cat("Saved: protein_filtered_80pct_outlier_removed.csv (",
    nrow(protein_filtered_reduced), "rows )\n")
cat("Saved: outlier_removed_sample_ids.rds\n")

