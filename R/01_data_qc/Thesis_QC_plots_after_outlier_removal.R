# Thesis-quality missingness and outlier QC plots
# Run the filtering script and outlier script first. This script uses their outputs.

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(pheatmap)

# ---- Set folders ----------------------------------------------------------
qc_dir <- "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Missing values and filtering NEW"
outlier_dir <- "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Missing values and filtering NEW"
plot_dir <- file.path(qc_dir, "thesis_qc_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

# The upstream scripts may have been run with either folder as the working directory.
find_input <- function(file_name) {
  candidates <- file.path(c(qc_dir, outlier_dir), file_name)
  found <- candidates[file.exists(candidates)]
  if (length(found) == 0) {
    stop("Cannot find ", file_name, ". Run the filtering and outlier scripts first, then place its output in either qc_dir or outlier_dir.")
  }
  found[1]
}

pre_file <- find_input("protein_filtered_80pct.csv")
post_file <- find_input("protein_filtered_80pct_outlier_removed.csv")
protein_file <- find_input("protein_column_names_post_filter.rds")
removed_file <- find_input("outlier_removed_sample_ids.rds")

protein_filtered <- read.csv(pre_file, check.names = FALSE)
protein_reduced <- read.csv(post_file, check.names = FALSE)
protein_cols <- readRDS(protein_file)
removed_ids <- readRDS(removed_file)

outlier_labels <- setNames(
  LETTERS[seq_along(removed_ids)],
  removed_ids
)

stopifnot(
  all(
    c(
      "SampleID",
      "case",
      "maternal_T1D",
      "corrected_visit_sample"
    ) %in% names(protein_reduced)
  )
)
stopifnot(all(protein_cols %in% names(protein_filtered)))
stopifnot(all(protein_cols %in% names(protein_reduced)))

# ---- Shared labels and requested colours ---------------------------------
case_colours <- c("Control" = "#4575b4", "Case" = "#d73027")
maternal_t1d_colours <- c("No maternal T1D" = "#91bfdb", "Maternal T1D" = "#fc8d59")
trimester_colours <- c(
  "T01" = "#ed693b",
  "T02" = "#43b097",
  "T03" = "#cb5299"
)

make_case_label <- function(x) {
  factor(ifelse(as.character(x) %in% c("1", "Case", "case"), "Case", "Control"), levels = c("Control", "Case"))
}

make_maternal_t1d_label <- function(x) {
  factor(ifelse(as.character(x) %in% c("1", "Yes", "yes", "TRUE", "True"), "Maternal T1D", "No maternal T1D"), levels = c("No maternal T1D", "Maternal T1D"))
}

make_trimester_label <- function(x) {
  factor(
    as.character(x),
    levels = c("T01", "T02", "T03")
  )
}

# ---- A. Final heatmap: after protein filtering AND outlier removal -------
final_mat <- protein_reduced %>% select(all_of(protein_cols)) %>% as.matrix()
rownames(final_mat) <- protein_reduced$SampleID
final_missing <- is.na(final_mat) * 1L
protein_order <- order(colMeans(final_missing), decreasing = TRUE)

row_annotation <- data.frame(
  Offspring_IA = make_case_label(protein_reduced$case),
  Maternal_T1D = make_maternal_t1d_label(protein_reduced$maternal_T1D),
  Trimester = make_trimester_label(
    protein_reduced$corrected_visit_sample
  ),
  row.names = protein_reduced$SampleID
)

annotation_colours <- list(
  Offspring_IA = case_colours,
  Maternal_T1D = maternal_t1d_colours,
  Trimester = trimester_colours
)

png(file.path(plot_dir, "Figure_QC_final_missingness_heatmap.png"), width = 3600, height = 2400, res = 300)
pheatmap(
  final_missing[, protein_order, drop = FALSE],
  color = c("#91bfdb", "#b2182b"),
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_rownames = FALSE, show_colnames = FALSE,
  annotation_row = row_annotation, annotation_colors = annotation_colours,
  legend_breaks = c(0, 1), legend_labels = c("Present", "Missing"),
  main = paste0("Missingness after protein filtering and outlier removal\n", nrow(final_mat), " retained samples x ", ncol(final_mat), " retained proteins"),
  fontsize = 10
)
dev.off()

# ---- B and C. Sample-level QC before removal -----------------------------
# This matches the intensity definition used in the existing outlier script.
pre_mat <- protein_filtered %>% select(all_of(protein_cols)) %>% as.matrix()
rownames(pre_mat) <- protein_filtered$SampleID

sample_qc <- protein_filtered %>%
  transmute(
    SampleID,
    Outlier_label = unname(outlier_labels[as.character(SampleID)]),
    Offspring_IA = make_case_label(case),
    Maternal_T1D = make_maternal_t1d_label(maternal_T1D),
    Removal_status = ifelse(SampleID %in% removed_ids, "Removed outlier", "Retained")
  ) %>%
  mutate(
    total_log2_intensity = rowSums(log2(pre_mat[SampleID, , drop = FALSE] + 1), na.rm = TRUE),
    pct_missing = rowMeans(is.na(pre_mat[SampleID, , drop = FALSE])) * 100
  ) %>%
  arrange(total_log2_intensity) %>%
  mutate(sample_rank = row_number())

intensity_cutoff <- mean(sample_qc$total_log2_intensity) - 3 * sd(sample_qc$total_log2_intensity)

p_intensity <- ggplot(sample_qc, aes(x = sample_rank, y = total_log2_intensity)) +
  geom_point(aes(fill = Offspring_IA, shape = Maternal_T1D), size = 3, colour = "#333333", stroke = 0.35) +
  geom_point(data = filter(sample_qc, Removal_status == "Removed outlier"), shape = 1, size = 4.7, colour = "#000000", stroke = 1) +
  geom_text_repel(data = filter(sample_qc, Removal_status == "Removed outlier"), aes(label = Outlier_label), size = 3, max.overlaps = Inf, min.segment.length = 0) +
  geom_hline(yintercept = intensity_cutoff, linetype = "dashed", colour = "#b2182b") +
  scale_fill_manual(values = case_colours, name = "Offspring IA outcome") +
  scale_shape_manual(values = c("No maternal T1D" = 22, "Maternal T1D" = 24), name = "Maternal T1D") +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        fill = unname(case_colours),
        colour = "#333333",
        size = 4
      )
    )
  ) +
  labs(title = "Sample total intensity before outlier removal", subtitle = "Black outlines and labels identify samples removed by the multi-method outlier rule", x = "Sample rank by total observed log2 intensity", y = "Total observed log2 intensity") +
  theme_classic(base_size = 12) + theme(legend.position = "top")

ggsave(file.path(plot_dir, "Figure_QC_total_intensity_outliers.png"), p_intensity, width = 10, height = 6, dpi = 600)

p_missingness_intensity <- ggplot(sample_qc, aes(x = total_log2_intensity, y = pct_missing)) +
  geom_point(aes(fill = Offspring_IA, shape = Maternal_T1D), size = 3, colour = "#333333", stroke = 0.35) +
  geom_point(data = filter(sample_qc, Removal_status == "Removed outlier"), shape = 1, size = 4.7, colour = "#000000", stroke = 1) +
  geom_text_repel(data = filter(sample_qc, Removal_status == "Removed outlier"), aes(label = Outlier_label), size = 3, max.overlaps = Inf, min.segment.length = 0) +
  scale_fill_manual(values = case_colours, name = "Offspring IA outcome") +
  scale_shape_manual(values = c("No maternal T1D" = 22, "Maternal T1D" = 24), name = "Maternal T1D") +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        fill = unname(case_colours),
        colour = "#333333",
        size = 4
      )
    )
  ) +
  labs(title = "Sample missingness and total intensity before outlier removal", subtitle = "Black outlines and labels identify samples removed by the multi-method outlier rule", x = "Total observed log2 intensity", y = "Missing values per sample (%)") +
  theme_classic(base_size = 12) + theme(legend.position = "top")

ggsave(file.path(plot_dir, "Figure_QC_missingness_vs_intensity_outliers.png"), p_missingness_intensity, width = 9, height = 6, dpi = 600)

# ---- D. Detection by case/control: documents the exact 80% rule ----------
detect_case <- colMeans(!is.na(pre_mat[protein_filtered$case == 1, , drop = FALSE])) * 100
detect_control <- colMeans(!is.na(pre_mat[protein_filtered$case == 0, , drop = FALSE])) * 100

detection_df <- data.frame(Protein = protein_cols, Case = detect_case, Control = detect_control) %>%
  mutate(Retained = ifelse(Case >= 80 | Control >= 80, "Retained", "Removed"))

p_detection <- ggplot(detection_df, aes(x = Control, y = Case, colour = Retained)) +
  geom_point(size = 2.1, alpha = 0.8) +
  geom_vline(xintercept = 80, linetype = "dashed", colour = "#555555") +
  geom_hline(yintercept = 80, linetype = "dashed", colour = "#555555") +
  scale_colour_manual(values = c("Retained" = "#2166ac", "Removed" = "#b2182b")) +
  coord_equal(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(title = "Protein detection by offspring IA outcome", subtitle = "Proteins were retained when detected in at least 80% of cases or controls", x = "Detected in controls (%)", y = "Detected in cases (%)", colour = "Filtering outcome") +
  theme_classic(base_size = 12) + theme(legend.position = "top")

ggsave(file.path(plot_dir, "Figure_QC_detection_by_outcome.png"), p_detection, width = 7, height = 6, dpi = 600)

write.csv(sample_qc, file.path(plot_dir, "sample_qc_summary_before_outlier_removal.csv"), row.names = FALSE)
message("Saved four thesis QC plots and a sample-QC summary to: ", plot_dir)