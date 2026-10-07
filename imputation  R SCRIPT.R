# Impute missing protein values and create diagnostic plots

library(tidyverse)
library(reshape2)
library(pcaMethods)
library(patchwork)
library(impute)
library(imputeLCMD)

set.seed(123)

# Load the filtered dataset and protein column names
protein_reduced <- read.csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Missing values and filtering NEW/protein_filtered_80pct_outlier_removed.csv",
  check.names = FALSE
)

prot_cols_expected <- readRDS(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Missing values and filtering NEW/protein_column_names_post_filter.rds"
)

missing_prot_cols <- setdiff(prot_cols_expected, names(protein_reduced))
if (length(missing_prot_cols) > 0) {
  stop("Protein columns are missing: ", paste(missing_prot_cols, collapse = ", "))
}

prot_cols <- prot_cols_expected
meta_cols <- setdiff(names(protein_reduced), prot_cols)
meta_qc <- protein_reduced %>% select(all_of(meta_cols))

prot_mat <- protein_reduced %>%
  select(all_of(prot_cols)) %>%
  as.matrix()
rownames(prot_mat) <- protein_reduced$SampleID

# Prepare log2 data for QRILC and kNN
log_mat <- log2(prot_mat + 1)
missing_mask <- is.na(log_mat)

cat("Samples:", nrow(log_mat), "| Proteins:", ncol(log_mat), "\n")
cat("Missing values before imputation:", sum(missing_mask), "\n")

# Plot colours
case_colours <- c("Control" = "#4575b4", "Case" = "#d73027")
plate_colours <- c("1" = "#6AA19E", "2" = "#F9F5E3", "3" = "#D97D8E")
rle_plate_colours <- c("1" = "#42B7B9", "2" = "#F2D14B", "3" = "#D691C1")
method_colours <- c(
  "Observed" = "#4d4d4d",
  "Min/2" = "#cc79a7",
  "kNN" = "#009e73",                         
  "QRILC" = "#0072b2"
)

# Create an RLE boxplot for each sample
plot_RLE <- function(mat, title, meta_df = meta_qc, color_by = "Plate.x") {
  protein_medians <- apply(mat, 2, median, na.rm = TRUE)
  rle_mat <- sweep(mat, 2, protein_medians, "-")

  rle_long <- reshape2::melt(
    as.matrix(rle_mat),
    varnames = c("SampleID", "Protein"),
    value.name = "RLE"
  ) %>%
    filter(!is.na(RLE)) %>%
    left_join(
      meta_df %>% select(SampleID, Group = all_of(color_by)),
      by = "SampleID"
    ) %>%
    mutate(Group = factor(Group))

  ggplot(rle_long, aes(x = SampleID, y = RLE, fill = Group)) +
    geom_boxplot(outlier.size = 0.15, outlier.alpha = 0.2, linewidth = 0.25) +
    geom_hline(yintercept = 0, colour = "#555555", linetype = "dashed") +
    scale_fill_manual(values = rle_plate_colours, name = color_by, drop = FALSE) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "top"
    ) +
    labs(title = title, x = "Samples", y = "Relative log2 abundance")
}

# Create a PCA plot; NIPALS is used before imputation to handle missing values
plot_PCA <- function(mat, title, meta_df = meta_qc) {
  pca_method <- if (anyNA(mat)) "nipals" else "svd"
  pca_res <- pca(
    as.matrix(mat),
    method = pca_method,
    nPcs = 2,
    scale = "uv",
    center = TRUE
  )

  pca_df <- as.data.frame(scores(pca_res)) %>%
    rownames_to_column("SampleID") %>%
    left_join(
      meta_df %>% select(SampleID, Plate.x, case),
      by = "SampleID"
    ) %>%
    mutate(
      Plate.x = factor(Plate.x),
      Case = factor(case, levels = c(0, 1), labels = c("Control", "Case"))
    )

  variance <- round(pca_res@R2[1:2] * 100, 1)

  ggplot(pca_df, aes(x = PC1, y = PC2, fill = Case, shape = Plate.x)) +
    geom_point(size = 2.8, alpha = 0.85, colour = "#333333", stroke = 0.3) +
    scale_fill_manual(values = case_colours, name = "Offspring IA") +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24), name = "Plate") +
    theme_classic(base_size = 11) +
    theme(legend.position = "top") +
    labs(
      title = title,
      x = paste0("PC1 (", variance[1], "%)"),
      y = paste0("PC2 (", variance[2], "%)")
    )
}

# Imputation method 1: half the minimum observed raw intensity per protein
impute_min_half <- function(mat) {
  for (j in seq_len(ncol(mat))) {
    missing <- is.na(mat[, j])
    if (any(missing)) {
      observed_min <- min(mat[, j], na.rm = TRUE)
      if (!is.finite(observed_min)) {
        stop("Protein column has no observed values: ", colnames(mat)[j])
      }
      mat[missing, j] <- observed_min / 2
    }
  }
  mat
}

imp_min <- log2(impute_min_half(prot_mat) + 1)

# Imputation methods 2 and 3 use the log2-transformed matrix
imp_knn <- t(impute.knn(t(log_mat), k = 10)$data)
imp_qrilc <- t(impute.QRILC(t(log_mat))[[1]])

imputed_list <- list(Min2 = imp_min, kNN = imp_knn, QRILC = imp_qrilc)

for (method in names(imputed_list)) {
  cat(method, "remaining missing values:", sum(is.na(imputed_list[[method]])), "\n")
}

# RLE plots
p_rle_before <- plot_RLE(log_mat, "Before imputation")
p_rle_min <- plot_RLE(imp_min, "After Min/2 imputation")
p_rle_knn <- plot_RLE(imp_knn, "After kNN imputation")
p_rle_qrilc <- plot_RLE(imp_qrilc, "After QRILC imputation")

ggsave("RLE_before_imputation.png", p_rle_before, width = 12, height = 5, dpi = 600)
ggsave("RLE_after_min2.png", p_rle_min, width = 12, height = 5, dpi = 600)
ggsave("RLE_after_kNN.png", p_rle_knn, width = 12, height = 5, dpi = 600)
ggsave("RLE_after_QRILC.png", p_rle_qrilc, width = 12, height = 5, dpi = 600)

p_rle_comparison <- (p_rle_before | p_rle_min) / (p_rle_knn | p_rle_qrilc) +
  plot_annotation(title = "Relative log abundance before and after imputation")
ggsave("RLE_comparison_panel.png", p_rle_comparison, width = 16, height = 10, dpi = 600)

# Compare imputed values with the observed values on the log2 scale
density_df <- bind_rows(
  tibble(value = as.numeric(imp_min[missing_mask]), method = "Min/2"),
  tibble(value = as.numeric(imp_knn[missing_mask]), method = "kNN"),
  tibble(value = as.numeric(imp_qrilc[missing_mask]), method = "QRILC"),
  tibble(value = as.numeric(log_mat[!missing_mask]), method = "Observed")
) %>%
  filter(is.finite(value)) %>%
  mutate(method = factor(method, levels = names(method_colours)))

p_density <- ggplot(density_df, aes(x = value, colour = method, fill = method)) +
  geom_density(alpha = 0.12, linewidth = 0.9, na.rm = TRUE) +
  scale_colour_manual(values = method_colours, name = "Values") +
  scale_fill_manual(values = method_colours, name = "Values") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top") +
  labs(
    title = "Observed and imputed protein abundance distributions",
    subtitle = "Imputed values are shown at positions that were missing before imputation",
    x = "log2(Intensity + 1)",
    y = "Density"
  )

ggsave("imputation_density_comparison.png", p_density, width = 10, height = 6, dpi = 600)

# PCA plots before and after imputation
p_pca_before <- plot_PCA(log_mat, "Before imputation")
p_pca_min <- plot_PCA(imp_min, "After Min/2 imputation")
p_pca_knn <- plot_PCA(imp_knn, "After kNN imputation")
p_pca_qrilc <- plot_PCA(imp_qrilc, "After QRILC imputation")

ggsave("PCA_before_imputation.png", p_pca_before, width = 8, height = 6, dpi = 600)
ggsave("PCA_after_min2.png", p_pca_min, width = 8, height = 6, dpi = 600)
ggsave("PCA_after_kNN.png", p_pca_knn, width = 8, height = 6, dpi = 600)
ggsave("PCA_after_QRILC.png", p_pca_qrilc, width = 8, height = 6, dpi = 600)

p_pca_comparison <- (p_pca_before | p_pca_min) / (p_pca_knn | p_pca_qrilc) +
  plot_annotation(title = "PCA comparison across imputation methods") +
  plot_layout(guides = "collect")
ggsave("PCA_imputation_comparison.png", p_pca_comparison, width = 14, height = 10, dpi = 600)

# Save imputed datasets and the protein column names
for (method in names(imputed_list)) {
  imputed_data <- meta_qc %>%
    left_join(
      as.data.frame(imputed_list[[method]]) %>% rownames_to_column("SampleID"),
      by = "SampleID"
    )
  write.csv(imputed_data, paste0("protein_imputed_", method, ".csv"), row.names = FALSE)
}

saveRDS(prot_cols, "protein_column_names_post_imputation.rds")
cat("Saved imputed datasets, QC plots, and protein_column_names_post_imputation.rds\n")
