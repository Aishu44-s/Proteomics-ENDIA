library(tidyverse)
library(reshape2)
library(limma)            # cyclic loess
library(vsn)               # VSN
library(preprocessCore)    # quantile
library(patchwork)
library(forcats)           # fct_reorder()

# Plot colours
case_colours <- c("Control" = "#4575b4", "Case" = "#d73027")
plate_colours <- c("1" = "#167C65", "2" = "#AFAC41", "3" = "#18608B")
rle_plate_colours <- c("1" = "#42B7B9", "2" = "#F2D14B", "3" = "#D691C1")
rle_visit_colours <- c("T01" = "#ed693b", "T02" = "#43b097", "T03" = "#cb5299")
norm_method_colours <- c(
  "Before (QRILC only)" = "#42B7B9",
  "VSN" = "#F2D14B",
  "Median" = "#D691C1",
  "Quantile" = "#D81B60",                              
  "CyclicLoess" = "#7B1FA2"
)

# 1. Load QRILC-imputed data (log2 scale) + metadata ------------------------

imputation_dir <- "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Imputation NEW"

protein_imputed <- read.csv(
  file.path(imputation_dir, "protein_imputed_QRILC.csv"),
  header = TRUE,
  check.names = FALSE
)

prot_cols_expected <- readRDS(
  file.path(imputation_dir, "protein_column_names_post_imputation.rds")
)

missing_prot_cols <- setdiff(prot_cols_expected, names(protein_imputed))
cat("Expected protein columns missing from protein_imputed (should be 0):",
    length(missing_prot_cols), "\n")
if (length(missing_prot_cols) > 0) print(missing_prot_cols)

prot_cols <- intersect(prot_cols_expected, names(protein_imputed))
meta_cols <- setdiff(names(protein_imputed), prot_cols)

meta_qc <- protein_imputed %>%
  dplyr::select(all_of(meta_cols)) %>%
  dplyr::rename(Plate = Plate.x) %>%
  dplyr::select(-any_of("Plate.y"))

imp_qrilc <- protein_imputed %>% dplyr::select(all_of(prot_cols)) %>% as.matrix()
rownames(imp_qrilc) <- protein_imputed$SampleID

cat("Samples:", nrow(imp_qrilc), " | Proteins:", ncol(imp_qrilc), "\n")
cat("NAs :", sum(is.na(imp_qrilc)), "\n")

# 2. Independent sample orders, one per grouping variable --------------------
# Plate: expect tightening after normalisation (batch correction working)
# Case / Timepoint: expect similar before/after (biology preserved)

order_by_var <- function(meta_df, var) {
  meta_df %>%
    mutate(SampleID = fct_reorder(SampleID, as.numeric(factor(.data[[var]])))) %>%
    arrange(SampleID) %>%
    pull(SampleID) %>%
    as.character()
}

order_plate <- order_by_var(meta_qc, "Plate")
order_case  <- order_by_var(meta_qc, "case")
order_visit <- order_by_var(meta_qc, "corrected_visit_sample")

# 3. Reusable plotting functions ---------------------------------------------
# plot_RLE orders and colours by the same group_var, by design

plot_RLE_panel <- function(
    mat,
    meta_df = meta_qc,
    group_var,
    sample_order
) {
  prot_medians <- apply(mat, 2, median, na.rm = TRUE)
  rle_mat <- sweep(mat, 2, prot_medians, "-")
  
  rle_df <- melt(
    as.matrix(rle_mat),
    varnames = c("SampleID", "Protein"),
    value.name = "RLE"
  ) %>%
    filter(!is.na(RLE)) %>%
    left_join(
      meta_df %>%
        dplyr::select(
          SampleID,
          case,
          corrected_visit_sample,
          Plate
        ),
      by = "SampleID"
    ) %>%
    mutate(
      SampleID = factor(SampleID, levels = sample_order),
      Group = case_when(
        group_var == "case" ~ factor(
          case,
          levels = c(0, 1),
          labels = c("Control", "Case")
        ),
        group_var == "Plate" ~ factor(Plate),
        group_var == "corrected_visit_sample" ~ factor(
          corrected_visit_sample,
          levels = c("T01", "T02", "T03")
        )
      )
    )
  
  legend_title <- switch(
    group_var,
    "case" = "Case/control",
    "Plate" = "Plate",
    "corrected_visit_sample" = "Visit"
  )
  
  x_label <- switch(
    group_var,
    "case" = "Samples ordered by case/control",
    "Plate" = "Samples ordered by plate",
    "corrected_visit_sample" = "Samples ordered by visit"
  )
  
  fill_colours <- switch(
    group_var,
    "case" = case_colours,
    "Plate" = rle_plate_colours,
    "corrected_visit_sample" = rle_visit_colours
  )
  
  ggplot(rle_df, aes(x = SampleID, y = RLE, fill = Group)) +
    geom_boxplot(
      outlier.size = 0.15,
      outlier.alpha = 0.2,
      linewidth = 0.25,
      width = 0.8
    ) +
    geom_hline(
      yintercept = 0,
      colour = "#555555",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_fill_manual(
      values = fill_colours,
      name = legend_title
    ) +
    labs(
      x = x_label,
      y = "Relative log2 expression"
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.title = element_text(size = 16),
      axis.text.y = element_text(size = 13),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "top",
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13)
    )
}
plot_density <- function(mat, title_str, meta_df = meta_qc) {
  dens_df <- melt(as.matrix(mat), varnames = c("SampleID", "Protein"),
                  value.name = "Intensity") %>%
    filter(!is.na(Intensity)) %>%
    left_join(meta_df %>% dplyr::select(SampleID, Plate), by = "SampleID") %>%
    mutate(Plate = factor(Plate))
  
  ggplot(dens_df, aes(x = Intensity, group = SampleID, colour = Plate)) +
    geom_density(linewidth = 0.35, alpha = 0.6) +
    scale_colour_manual(values = plate_colours, name = "Plate") +
    theme_classic(base_size = 12) +
    theme(panel.background = element_rect(fill = "#EEEEEE", colour = NA)) +
    theme(legend.position = "top") +
    labs(title = title_str, subtitle = "Each line = one sample | Should overlap tightly after normalisation",
         x = "log2 intensity", y = "Density")
}

plot_PCA <- function(mat, title_str, meta_df = meta_qc) {
  pca_res <- prcomp(mat, center = TRUE, scale. = FALSE)
  pca_df <- as.data.frame(pca_res$x[, 1:2]) %>%
    rownames_to_column("SampleID") %>%
    left_join(meta_df %>% dplyr::select(SampleID, Plate, case), by = "SampleID") %>%
    mutate(Plate = factor(Plate),
           Case  = factor(case, levels = c(0, 1), labels = c("Control", "Case")))
  
  var1 <- round(summary(pca_res)$importance[2, 1] * 100, 1)
  var2 <- round(summary(pca_res)$importance[2, 2] * 100, 1)
  
  ggplot(pca_df, aes(x = PC1, y = PC2, fill = Case, shape = Plate)) +
    geom_point(size = 2.8, alpha = 0.9, colour = "#333333", stroke = 0.3) +
    stat_ellipse(aes(group = Plate), type = "norm", linetype = "dashed",
                 linewidth = 0.5, alpha = 0.7) +
    scale_fill_manual(values = case_colours, name = "IA Status") +
    scale_shape_manual(values = c("1" = 21, "2" = 22, "3" = 24), name = "Plate") +
    theme_classic(base_size = 12) +
    labs(title = title_str, x = paste0("PC1 (", var1, "%)"), y = paste0("PC2 (", var2, "%)"))
}

# 4. Run normalisation methods -----------------------------------------------

raw_scale <- 2^imp_qrilc - 1
raw_scale[raw_scale < 0] <- 0

vsn_fit  <- vsn2(t(as.matrix(raw_scale)))
norm_vsn <- as.data.frame(t(predict(vsn_fit, t(as.matrix(raw_scale)))))
rownames(norm_vsn) <- rownames(imp_qrilc)
colnames(norm_vsn)  <- colnames(imp_qrilc)

sample_medians <- apply(imp_qrilc, 1, median, na.rm = TRUE)
norm_median <- sweep(imp_qrilc, 1, sample_medians, "-")

norm_quant <- as.data.frame(t(normalize.quantiles(t(as.matrix(imp_qrilc)))))
rownames(norm_quant) <- rownames(imp_qrilc)
colnames(norm_quant)  <- colnames(imp_qrilc)

norm_cyclic <- as.data.frame(t(normalizeCyclicLoess(t(as.matrix(imp_qrilc)), method = "fast")))
rownames(norm_cyclic) <- rownames(imp_qrilc)
colnames(norm_cyclic)  <- colnames(imp_qrilc)

norm_methods <- list("VSN" = norm_vsn, "Median" = norm_median,
                     "Quantile" = norm_quant, "CyclicLoess" = norm_cyclic)

cat("Normalisation complete for:", paste(names(norm_methods), collapse = ", "), "\n")

# 5. Combined RLE plot: before and after VSN normalisation -------------------

# Before normalisation
p_before_case <- plot_RLE_panel(
  imp_qrilc,
  group_var = "case",
  sample_order = order_case
) +
  labs(title = "Before normalisation: case/control")

p_before_plate <- plot_RLE_panel(
  imp_qrilc,
  group_var = "Plate",
  sample_order = order_plate
) +
  labs(title = "Before normalisation: plate")

p_before_visit <- plot_RLE_panel(
  imp_qrilc,
  group_var = "corrected_visit_sample",
  sample_order = order_visit
) +
  labs(title = "Before normalisation: visit")

# After VSN normalisation
p_vsn_case <- plot_RLE_panel(
  norm_vsn,
  group_var = "case",
  sample_order = order_case
) +
  labs(title = "VSN normalisation: case/control")

p_vsn_plate <- plot_RLE_panel(
  norm_vsn,
  group_var = "Plate",
  sample_order = order_plate
) +
  labs(title = "VSN normalisation: plate")

p_vsn_visit <- plot_RLE_panel(
  norm_vsn,
  group_var = "corrected_visit_sample",
  sample_order = order_visit
) +
  labs(title = "VSN normalisation: visit")

p_rle_before_after <- (
  (p_before_case  | p_vsn_case) /
    (p_before_plate | p_vsn_plate) /
    (p_before_visit | p_vsn_visit)
) +
  plot_layout(
    guides = "collect",
    heights = c(1, 1, 1)
  ) +
  plot_annotation(
    title = "RLE before and after VSN normalisation",
    subtitle = "Samples ordered by case/control, plate and visit",
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 28,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 18,
        hjust = 0.5
      )
    )
  ) &
  theme(
    legend.position = "top",
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14),
    plot.title = element_text(
      face = "bold",
      size = 18,
      hjust = 0.5
    )
  )

print(p_rle_before_after)

vsn_plot_dir <- paste0(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/",
  "NEW ANALYSIS/Normalisation NEW/VSN PLOTS"
)

ggsave(
  filename = file.path(
    vsn_plot_dir,
    "RLE_before_after_VSN_case_plate_visit.png"
  ),
  plot = p_rle_before_after,
  width = 18,
  height = 24,
  dpi = 400,
  bg = "white"
)

# Density plot before normalisation ------------------------------------------

p_density_before <- plot_density(
  imp_qrilc,
  "Density - Before normalisation (QRILC-imputed)"
)

ggsave(
  filename = file.path(
    vsn_plot_dir,
    "norm_density_before_normalisation.png"
  ),
  plot = p_density_before,
  width = 10,
  height = 6,
  dpi = 600,
  bg = "white"
)

# Keep density and PCA figures for every normalisation method
for (method_name in names(norm_methods)) {
  
  mat <- as.matrix(norm_methods[[method_name]])
  
  p_dens <- plot_density(
    mat,
    paste("Density -", method_name, "Normalisation")
  )
  
  ggsave(
    paste0("norm_density_", method_name, ".png"),
    plot = p_dens,
    width = 10,
    height = 6,
    dpi = 600
  )
  
  p_pca <- plot_PCA(
    mat,
    paste("PCA -", method_name, "Normalisation")
  )
  
  ggsave(
    paste0("norm_PCA_", method_name, ".png"),
    plot = p_pca,
    width = 8,
    height = 6,
    dpi = 600
  )
}

# 6. PMAD comparison (lower = better) -----------------------------------------

compute_PMAD <- function(mat, label) {
  data.frame(PMAD = apply(mat, 2, function(x) mad(x, na.rm = TRUE)), Method = label)
}

pmad_all <- bind_rows(
  compute_PMAD(imp_qrilc,     "Before (QRILC only)"),
  compute_PMAD(norm_vsn,      "VSN"),
  compute_PMAD(norm_median,   "Median"),
  compute_PMAD(norm_quant,    "Quantile"),
  compute_PMAD(norm_cyclic,   "CyclicLoess")
)

pmad_all$Method <- factor(pmad_all$Method,
                          levels = c("Before (QRILC only)", "VSN", "Median", "Quantile", "CyclicLoess"))

pmad_summary <- pmad_all %>%
  group_by(Method) %>%
  summarise(median_PMAD = round(median(PMAD), 4), .groups = "drop") %>%
  arrange(median_PMAD)

cat("\n=== PMAD SUMMARY (lower = better) ===\n")
print(pmad_summary)

p_pmad <- ggplot(pmad_all, aes(x = Method, y = PMAD, fill = Method)) +
  geom_boxplot(outlier.size = 0.3, outlier.alpha = 0.3) +
  geom_point(data = pmad_summary, aes(x = Method, y = median_PMAD),
             colour = "black", size = 3.5, shape = 18) +
  scale_fill_manual(values = norm_method_colours) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1)) +
  labs(title = "PMAD - Normalisation Method Comparison",
       subtitle = "Lower = better | Diamond = median PMAD per method",
       x = NULL, y = "Median Absolute Deviation")

ggsave("norm_PMAD_comparison.png", plot = p_pmad, width = 10, height = 6, dpi = 600)
cat("Saved: norm_PMAD_comparison.png\n")

# 7. Save normalised matrices --------------------------------------------------

for (method_name in names(norm_methods)) {
  out <- meta_qc %>%
    left_join(as.data.frame(norm_methods[[method_name]]) %>% rownames_to_column("SampleID"),
              by = "SampleID")
  write.csv(out, paste0("protein_normalised_", method_name, ".csv"), row.names = FALSE)
}
cat("Saved: protein_normalised_[VSN/Median/Quantile/CyclicLoess].csv\n")

saveRDS(prot_cols, "protein_column_names_post_normalisation.rds")
cat("Saved: protein_column_names_post_normalisation.rds\n")
