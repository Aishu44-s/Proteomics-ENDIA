library(tidyverse)
library(reshape2)
library(pcaMethods)
library(patchwork)
library(impute)       # kNN
library(imputeLCMD)   # QRILC

set.seed(123)   # QRILC has a random component — keep this for reproducibility

# 1. Load data, log2-transform ------------------------------------------

protein_reduced <- read.csv("protein_filtered_80pct_outlier_removed.csv", header = TRUE)

prot_cols_expected <- readRDS("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/protein_column_names_post_filter.rds")

missing_prot_cols <- setdiff(prot_cols_expected, names(protein_reduced))
cat("Expected protein columns missing from protein_reduced (should be 0):",
    length(missing_prot_cols), "\n")
if (length(missing_prot_cols) > 0) print(missing_prot_cols)

prot_cols <- intersect(prot_cols_expected, names(protein_reduced))
meta_cols <- setdiff(names(protein_reduced), prot_cols)

meta_qc  <- protein_reduced %>% dplyr::select(all_of(meta_cols))
prot_mat <- protein_reduced %>% dplyr::select(all_of(prot_cols)) %>% as.matrix()
rownames(prot_mat) <- protein_reduced$SampleID

log_mat <- log2(prot_mat + 1)   # +1 avoids log2(0) = -Inf

cat("Samples:", nrow(log_mat), " | Proteins:", ncol(log_mat), "\n")
cat("NAs before imputation:", sum(is.na(log_mat)), "\n")

# 2. Reusable plotting functions -----------------------------------------

# RLE: flags samples/plates whose overall level is shifted away from 0
plot_RLE <- function(mat, title, meta_df = meta_qc, color_by = "Plate.x") {
  prot_medians <- apply(mat, 2, median, na.rm = TRUE)
  rle_mat <- sweep(mat, 2, prot_medians, "-")
  
  rle_long <- melt(as.matrix(rle_mat), varnames = c("SampleID", "Protein"),
                   value.name = "RLE") %>%
    filter(!is.na(RLE)) %>%
    left_join(meta_df %>% dplyr::select(SampleID, Group = all_of(color_by)),
              by = "SampleID") %>%
    mutate(Group = factor(Group))
  
  ggplot(rle_long, aes(x = SampleID, y = RLE, fill = Group)) +
    geom_boxplot(outlier.size = 0.15, outlier.alpha = 0.2, linewidth = 0.25) +
    geom_hline(yintercept = 0, colour = "red", linetype = "dashed", linewidth = 0.7) +
    scale_fill_brewer(palette = "Set2", name = color_by) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "top") +
    labs(title = title, x = "Samples", y = "Relative log2 expression")
}

plot_PCA <- function(mat, title, meta_df = meta_qc) {
  pca_res <- pca(as.matrix(mat), method = "svd", nPcs = 2, scale = "uv", center = TRUE)
  
  pca_df <- as.data.frame(scores(pca_res)) %>%
    rownames_to_column("SampleID") %>%
    left_join(meta_df %>% dplyr::select(SampleID, Plate.x, case), by = "SampleID") %>%
    mutate(Plate.x = factor(Plate.x),
           Case  = factor(case, levels = c(0, 1), labels = c("Control", "Case")))
  
  var1 <- round(pca_res@R2[1] * 100, 1)
  var2 <- round(pca_res@R2[2] * 100, 1)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, colour = Plate.x, shape = Case)) +
    geom_point(size = 2.5, alpha = 0.85) +
    scale_colour_brewer(palette = "Set2") +
    theme_classic(base_size = 11) +
    labs(title = title, x = paste0("PC1 (", var1, "%)"), y = paste0("PC2 (", var2, "%)"))
}

# 3. Run the three imputation methods -------------------------------------

missing_mask <- is.na(log_mat)

# Method 1: min/2 per protein
impute_min_half <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    mat[is.na(mat[, j]), j] <- min(mat[, j], na.rm = TRUE) / 2
  }
  mat
}
imp_min <- impute_min_half(log_mat)

# Method 2: kNN (k = 10), needs proteins as rows
imp_knn <- t(impute.knn(t(log_mat), k = 10)$data)

# Method 3: QRILC, needs proteins as rows
imp_qrilc <- t(impute.QRILC(t(log_mat))[[1]])

imputed_list <- list("Min2" = imp_min, "kNN" = imp_knn, "QRILC" = imp_qrilc)

for (m in names(imputed_list)) {
  cat("Method:", m, "| NAs remaining:", sum(is.na(imputed_list[[m]])), "\n")
}

# 4. RLE plots — before AND after each method ------------------------------

p_rle_before <- plot_RLE(log_mat,   "RLE — Before Imputation")
ggsave("RLE_before_imputation.png", plot = p_rle_before, width = 12, height = 5, dpi = 600)

p_rle_min <- plot_RLE(imp_min,      "RLE — After Min/2 per Protein")
ggsave("RLE_after_min2.png", plot = p_rle_min, width = 12, height = 5, dpi = 600)

p_rle_knn <- plot_RLE(imp_knn,      "RLE — After kNN")
ggsave("RLE_after_kNN.png", plot = p_rle_knn, width = 12, height = 5, dpi = 600)

p_rle_qrilc <- plot_RLE(imp_qrilc,  "RLE — After QRILC")
ggsave("RLE_after_QRILC.png", plot = p_rle_qrilc, width = 12, height = 5, dpi = 600)

p_rle_combined <- (p_rle_before | p_rle_min) / (p_rle_knn | p_rle_qrilc) +
  plot_annotation(title = "RLE Comparison: Before vs. After Each Imputation Method")
ggsave("RLE_comparison_panel.png", plot = p_rle_combined, width = 16, height = 10, dpi = 600)

cat("Saved RLE plots (individual + combined panel)\n")

# 5. Density plot — where do imputed values land, per method? ---------------

density_df <- bind_rows(
  data.frame(value = as.numeric(as.matrix(imp_min)[missing_mask]),   method = "Min/2 per protein"),
  data.frame(value = as.numeric(as.matrix(imp_knn)[missing_mask]),   method = "kNN"),
  data.frame(value = as.numeric(as.matrix(imp_qrilc)[missing_mask]), method = "QRILC"),
  data.frame(value = as.numeric(log_mat[!missing_mask]),             method = "Observed (reference)")
)

p_density <- ggplot(density_df, aes(x = value, colour = method, fill = method)) +
  geom_density(alpha = 0.15, linewidth = 0.9) +
  scale_colour_manual(values = c("Observed (reference)" = "black", "Min/2 per protein" = "#E41A1C",
                                 "kNN" = "#4DAF4A", "QRILC" = "#377EB8")) +
  scale_fill_manual(values = c("Observed (reference)" = "grey80", "Min/2 per protein" = "#E41A1C",
                               "kNN" = "#4DAF4A", "QRILC" = "#377EB8")) +
  theme_classic(base_size = 13) +
  labs(title = "Where Do Imputed Values Land? — Method Comparison",
       subtitle = "A good MNAR method places imputed values LEFT of the observed distribution",
       x = "log2(Intensity + 1)", y = "Density", colour = "Method", fill = "Method")

ggsave("imputation_density_comparison.png", plot = p_density, width = 10, height = 6, dpi = 600)
cat("Saved: imputation_density_comparison.png\n")

# 6. PCA plots — individual + combined comparison panel ---------------------

p_pca_before <- plot_PCA(log_mat,   "PCA — Before Imputation")
ggsave("PCA_before_imputation.png", plot = p_pca_before, width = 8, height = 6, dpi = 600)

p_pca_min <- plot_PCA(imp_min,      "PCA — After Min/2")
ggsave("PCA_after_min2.png", plot = p_pca_min, width = 8, height = 6, dpi = 600)

p_pca_knn <- plot_PCA(imp_knn,      "PCA — After kNN")
ggsave("PCA_after_kNN.png", plot = p_pca_knn, width = 8, height = 6, dpi = 600)

p_pca_qrilc <- plot_PCA(imp_qrilc,  "PCA — After QRILC")
ggsave("PCA_after_QRILC.png", plot = p_pca_qrilc, width = 8, height = 6, dpi = 600)

p_pca_combined <- (p_pca_before | p_pca_min) / (p_pca_knn | p_pca_qrilc) +
  plot_annotation(title = "PCA Comparison Across Imputation Methods") +
  plot_layout(guides = "collect")
ggsave("PCA_imputation_comparison.png", plot = p_pca_combined, width = 14, height = 10, dpi = 600)

cat("Saved PCA plots (individual + combined panel)\n")

# 7. Save each imputed matrix ------------------------------------------------

for (m in names(imputed_list)) {
  out <- meta_qc %>%
    left_join(as.data.frame(imputed_list[[m]]) %>% rownames_to_column("SampleID"),
              by = "SampleID")
  write.csv(out, paste0("protein_imputed_", m, ".csv"), row.names = FALSE)
}
cat("Saved: protein_imputed_[Min2/kNN/QRILC].csv\n")

saveRDS(prot_cols, "protein_column_names_post_imputation.rds")
cat("Saved: protein_column_names_post_imputation.rds\n")