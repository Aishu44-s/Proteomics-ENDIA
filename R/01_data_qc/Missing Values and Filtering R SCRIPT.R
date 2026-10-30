library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(stringr)

# 1. Load the data ---------------------------------------------------------

#merged <- read.csv("merged_meta_protein.csv", header = TRUE)
#unique_samples <- merged %>%
#  distinct(SampleID.2, .keep_all = TRUE)

cat("Unique samples used for QC:", nrow(unique_samples), "\n")

protein_col_names <- readRDS("protein_column_names.rds")

missing_prot_cols <- setdiff(protein_col_names, names(unique_samples))
cat("Expected protein columns missing from unique_samples (should be 0):",
    length(missing_prot_cols), "\n")
if (length(missing_prot_cols) > 0) print(missing_prot_cols)

prot_cols <- intersect(protein_col_names, names(unique_samples))
cat("Number of protein columns detected:", length(prot_cols), "\n")

intensity_mat <- unique_samples %>%
  select(all_of(prot_cols)) %>%
  as.matrix()
rownames(intensity_mat) <- unique_samples$SampleID

# 2. Protein intensity overview -----------------------------------------

intensity_long <- as.data.frame(intensity_mat) %>%
  rownames_to_column("SampleID") %>%
  pivot_longer(-SampleID, names_to = "Protein", values_to = "Intensity")

ggplot(intensity_long, aes(x = SampleID, y = log2(Intensity + 1))) +
  geom_boxplot(outlier.size = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  labs(title = "Log2 protein intensity distribution per sample",
       x = "Sample", y = "log2(Intensity + 1)")

ggsave("intensity_boxplot_per_sample.png", width = 12, height = 5)

ggplot(intensity_long, aes(x = log2(Intensity + 1))) +
  geom_density(fill = "#E94E77", alpha = 0.4) +
  theme_minimal() +
  labs(title = "Overall log2 protein intensity distribution",
       x = "log2(Intensity + 1)", y = "Density")

ggsave("intensity_density_overall.png", width = 7, height = 5)

# Protein abundance rank plot
protein_abundance <- data.frame(
  Protein             = prot_cols,
  mean_log2_intensity = log2(colMeans(intensity_mat, na.rm = TRUE) + 1)
) %>%
  arrange(desc(mean_log2_intensity)) %>%
  mutate(rank = row_number())

ggplot(protein_abundance, aes(x = rank, y = mean_log2_intensity)) +
  geom_point(size = 0.8, alpha = 0.6, colour = "#E94E77") +
  theme_minimal() +
  labs(title = "Protein abundance rank plot",
       subtitle = paste0("n = ", nrow(protein_abundance), " proteins, ranked by mean log2 intensity"),
       x = "Protein rank (highest to lowest abundance)",
       y = "Mean log2(Intensity + 1) (observed values only)")

ggsave("protein_abundance_rank_plot.png", width = 8, height = 5, dpi = 600)

# Total intensity per sample
plot_total_intensity <- function(color_var, colours, legend_title = color_var) {
  total_intensity <- data.frame(
    SampleID = unique_samples$SampleID,
    total = rowSums(intensity_mat, na.rm = TRUE),
    Group = factor(unique_samples[[color_var]])
  ) %>%
    arrange(total) %>%
    mutate(SampleID = factor(SampleID, levels = SampleID))
  
  ggplot(total_intensity, aes(x = SampleID, y = total, fill = Group)) +
    geom_col(width = 1) +
    scale_fill_manual(values = colours, name = legend_title) +
    theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(
      title = paste("Total intensity per sample by", legend_title),
      x = "Sample (ranked by total intensity)",
      y = "Total intensity"
    )
}

p_intensity_case <- plot_total_intensity(
  "case",
  colours = c("0" = "#4575b4", "1" = "#d73027"),
  legend_title = "Offspring IA outcome"
)

p_intensity_plate <- plot_total_intensity(
  "Plate.x",
  colours = c("1" = "#5F8A8A", "2" = "#B86B7B", "3" = "#8FA374"),
  legend_title = "Plate"
)

p_intensity_trimester <- plot_total_intensity(
  "corrected_visit_sample",
  colours = c("T01" = "#CCBB44", "T02" = "#EE6677", "T03" = "#66CCEE"),
  legend_title = "Trimester"
)

ggsave("total_intensity_by_case.png", p_intensity_case, width = 10, height = 6, dpi = 600)
ggsave("total_intensity_by_plate.png", p_intensity_plate, width = 10, height = 6, dpi = 600)
ggsave("total_intensity_by_trimester.png", p_intensity_trimester, width = 10, height = 6, dpi = 600)

# 3. Missingness per protein ---------------------------------------------

missing_summary <- data.frame(
  Protein     = prot_cols,
  n_missing   = colSums(is.na(intensity_mat)),
  pct_missing = colMeans(is.na(intensity_mat)) * 100
) %>%
  arrange(desc(pct_missing))

cat("\nMissingness summary (top 10 proteins by % missing):\n")
print(head(missing_summary, 10))

ggplot(missing_summary, aes(x = pct_missing)) +
  geom_histogram(binwidth = 5, fill = "steelblue", colour = "white") +
  geom_vline(xintercept = 80, colour = "red", linetype = "dashed", linewidth = 1) +
  theme_minimal() +
  labs(title = "Distribution of % missing values across proteins",
       subtitle = "Red line = 80% missing cutoff",
       x = "% missing", y = "Number of proteins")

ggsave("missingness_histogram.png", width = 7, height = 5)

# 4. Missingness heatmap ---------------------------------------------------

missing_mat <- ifelse(is.na(intensity_mat), 1, 0)
missing_mat <- missing_mat[, order(-missing_summary$pct_missing[match(colnames(missing_mat), missing_summary$Protein)])]

pheatmap(
  missing_mat,
  color = c("#81E6D9", "#2D3748"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  legend_breaks = c(0, 1),
  legend_labels = c("Present", "Missing"),
  main = "Missingness pattern: samples (rows) x proteins (columns)",
  filename = "missingness_heatmap.png",
  width = 10, height = 7, dpi  = 600
)

# 5. Filter: keep proteins present in >=80% of samples in EITHER group -----

case_ids <- unique_samples$SampleID[unique_samples$case == 1]
ctrl_ids <- unique_samples$SampleID[unique_samples$case == 0]

prot_case <- intensity_mat[rownames(intensity_mat) %in% case_ids, , drop = FALSE]
prot_ctrl <- intensity_mat[rownames(intensity_mat) %in% ctrl_ids, , drop = FALSE]

valid_case <- colMeans(!is.na(prot_case))
valid_ctrl <- colMeans(!is.na(prot_ctrl))

keep <- valid_case >= 0.8 | valid_ctrl >= 0.8

proteins_to_keep <- prot_cols[keep]
proteins_removed <- prot_cols[!keep]

cat("\nProteins before filtering:", length(prot_cols), "\n")
cat("Proteins removed (<80% present in both groups):", length(proteins_removed), "\n")
cat("Proteins kept (>=80% present in at least one group):", length(proteins_to_keep), "\n")

meta_cols <- setdiff(names(unique_samples), prot_cols)

protein_filtered <- unique_samples %>%
  select(all_of(meta_cols), all_of(proteins_to_keep))

write.csv(protein_filtered, "protein_filtered_80pct.csv", row.names = FALSE)
cat("\nSaved: protein_filtered_80pct.csv (",
    nrow(protein_filtered), "rows x", ncol(protein_filtered), "cols )\n")

saveRDS(proteins_to_keep, "protein_column_names_post_filter.rds")
cat("Saved: protein_column_names_post_filter.rds\n")
# 6. Missingness heatmap — after filtering ---------------------------------

intensity_mat_filtered <- intensity_mat[, proteins_to_keep, drop = FALSE]

missing_mat_filtered <- ifelse(is.na(intensity_mat_filtered), 1, 0)

pheatmap(
  missing_mat_filtered,
  color = c("#81E6D9", "#2D3748"),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  legend_breaks = c(0, 1),
  legend_labels = c("Present", "Missing"),
  main = "Missingness pattern after filtering: samples (rows) x proteins (columns)",
  filename = "missingness_heatmap_post_filter.png",
  width = 10, height = 7, dpi  = 600
)
cat("Saved: missingness_heatmap_post_filter.png\n")