#===================================================================
# Missingness Pattern Diagnostics: Case/Control & Trimester
#
# *** THIS SCRIPT DOES NOT FILTER ANY PROTEINS ***
# The 80%-overall-missingness filter (Missing_Values_and_Filtering_
# R_SCRIPT_FIXED.R) already decided which proteins to keep — that
# stays your one canonical filtering step, and matches your primer
# (0 proteins removed).
#
# This script answers a DIFFERENT question on the proteins you
# already kept: does missingness pattern systematically differ by
# case/control status or by trimester? That matters for interpreting
# your MAR/MNAR conclusions and for flagging any protein whose
# apparent "signal" might actually just be a missingness artifact.
#
# Uses the same canonical protein list + files as the rest of your
# fixed pipeline — NOT the raw files re-read from scratch, and NOT
# any group-wise filter.
#===================================================================

library(tidyverse)
library(ggrepel)
library(pheatmap)

#-------------------------------------------------------------------
# 1. Load data — same canonical source as the rest of the pipeline
#-------------------------------------------------------------------

#protein_filtered <- read.csv("protein_filtered_80pct.csv", header = TRUE)
#prot_cols <- readRDS("protein_column_names_post_filter.rds")

#meta_clean <- protein_filtered %>%
  select(-all_of(prot_cols)) %>%
  rename(Plate = Plate.x)

prot_mat <- protein_filtered %>%
  select(all_of(prot_cols)) %>%
  as.matrix()
rownames(prot_mat) <- protein_filtered$SampleID

cat("Samples:", nrow(prot_mat), "| Proteins:", ncol(prot_mat), "\n")
cat("(This should be 469 proteins, 0 removed by the overall filter)\n")

#-------------------------------------------------------------------
# 2. Per-protein missingness: case vs control
#-------------------------------------------------------------------

case_ids <- meta_clean$SampleID[meta_clean$case == 1]
ctrl_ids <- meta_clean$SampleID[meta_clean$case == 0]

miss_case <- colMeans(is.na(prot_mat[rownames(prot_mat) %in% case_ids, ])) * 100
miss_ctrl <- colMeans(is.na(prot_mat[rownames(prot_mat) %in% ctrl_ids, ])) * 100

group_miss_df <- data.frame(
  Protein      = prot_cols,
  pct_miss_case = miss_case,
  pct_miss_ctrl = miss_ctrl
) %>%
  mutate(
    abs_diff   = abs(pct_miss_case - pct_miss_ctrl),
    short_name = gsub(".*\\.", "", Protein)
  )

# Flag proteins with a notable case/control missingness gap —
# DIAGNOSTIC LABEL ONLY, does not remove anything
flagged_group_diff <- group_miss_df %>%
  filter(abs_diff > 15, pmax(pct_miss_case, pct_miss_ctrl) > 5) %>%
  arrange(desc(abs_diff))

cat("\nProteins with >15pp missingness difference between case/control:",
    nrow(flagged_group_diff), "\n")
if (nrow(flagged_group_diff) > 0) print(flagged_group_diff)

p_group_scatter <- ggplot(group_miss_df, aes(x = pct_miss_ctrl, y = pct_miss_case)) +
  geom_point(alpha = 0.5, colour = "steelblue", size = 1.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_text_repel(data = flagged_group_diff, aes(label = short_name),
                  size = 2.5, max.overlaps = 20) +
  labs(title = "% Missing per Protein: Case vs Control",
       subtitle = paste0(nrow(flagged_group_diff),
                         " proteins with >15pp difference (labelled) — flagged for awareness, not removed"),
       x = "% Missing in Controls", y = "% Missing in Cases") +
  theme_minimal(base_size = 12)

ggsave("missingness_case_vs_control.png", p_group_scatter, width = 8, height = 7, dpi = 600)
write.csv(flagged_group_diff, "flagged_case_control_missingness_diff.csv", row.names = FALSE)

#-------------------------------------------------------------------
# 3. Missingness heatmap — samples annotated by case, Plate, trimester
#-------------------------------------------------------------------

miss_binary <- is.na(prot_mat) * 1

annot_row <- data.frame(
  Case      = factor(meta_clean$case, labels = c("Control", "Case")),
  Plate     = factor(meta_clean$Plate),
  Trimester = factor(meta_clean$corrected_visit_sample),
  row.names = meta_clean$SampleID
)

annot_colors <- list(
  Case      = c(Control = "#4393C3", Case = "#D6604D"),
  Plate     = c("1" = "#66C2A5", "2" = "#FC8D62", "3" = "#8DA0CB"),
  Trimester = c(T01 = "#FEE090", T02 = "#ABD9E9", T03 = "#4575B4")
)

# Top 100 most-missing proteins, for readability
top_miss_idx <- order(colMeans(is.na(prot_mat)), decreasing = TRUE)[1:min(100, ncol(prot_mat))]

pdf("missingness_heatmap_annotated.pdf", width = 14, height = 10)
pheatmap(
  miss_binary[, top_miss_idx],
  color             = c("white", "#D6604D"),
  cluster_rows      = TRUE,
  cluster_cols      = TRUE,
  show_rownames     = FALSE,
  show_colnames     = FALSE,
  annotation_row    = annot_row,
  annotation_colors = annot_colors,
  main = "Missingness Pattern — Top 100 Most-Missing Proteins\nAnnotated by Case, Plate, Trimester",
  fontsize = 8
)
dev.off()
cat("Saved: missingness_heatmap_annotated.pdf\n")

#-------------------------------------------------------------------
# 4. Trimester-level missingness — density + Kruskal-Wallis
#-------------------------------------------------------------------

trimester_miss_long <- lapply(sort(unique(meta_clean$corrected_visit_sample)), function(v) {
  idx <- which(meta_clean$corrected_visit_sample == v)
  data.frame(
    Trimester = v,
    pct_miss  = colMeans(is.na(prot_mat[idx, , drop = FALSE])) * 100
  )
}) %>% bind_rows()

p_trimester <- ggplot(trimester_miss_long, aes(x = pct_miss, fill = Trimester)) +
  geom_density(alpha = 0.45, linewidth = 0.7) +
  labs(title = "Trimester-Level Missingness Distribution (per protein)",
       subtitle = "Separation between trimesters = a visit-related (MAR) component",
       x = "% Missing per Protein", y = "Density") +
  theme_classic(base_size = 13) +
  theme(legend.position = "top")

ggsave("missingness_by_trimester_density.png", p_trimester, width = 8, height = 5, dpi = 150)

# Does per-SAMPLE missingness differ across trimesters? (matches the
# same logic your MAR logistic regression already tests, shown here
# as a simple non-parametric cross-check)
sample_miss_pct <- rowMeans(is.na(prot_mat)) * 100
trimester_kw <- kruskal.test(sample_miss_pct ~ factor(meta_clean$corrected_visit_sample))

cat("\n--- Trimester-Level Missingness (Kruskal-Wallis) ---\n")
cat("p =", signif(trimester_kw$p.value, 4), "\n")
if (trimester_kw$p.value < 0.05) {
  cat("RESULT: Missingness differs significantly across trimesters (matches your primer's finding for T03)\n")
} else {
  cat("RESULT: No significant trimester-level difference in missingness\n")
}

# Same check, by Plate — cross-checks your primer's Plate 2 finding
plate_kw <- kruskal.test(sample_miss_pct ~ factor(meta_clean$Plate))
cat("\n--- Plate-Level Missingness (Kruskal-Wallis) ---\n")
cat("p =", signif(plate_kw$p.value, 4), "\n")

#-------------------------------------------------------------------
# 5. Case-level missingness — cross-checks your primer's finding
#    that case status was NOT a significant MAR driver
#-------------------------------------------------------------------

case_kw <- kruskal.test(sample_miss_pct ~ factor(meta_clean$case))
cat("\n--- Case/Control-Level Missingness (Kruskal-Wallis) ---\n")
cat("p =", signif(case_kw$p.value, 4), "\n")
if (case_kw$p.value >= 0.05) {
  cat("RESULT: No significant difference — consistent with your primer (case was NOT a MAR driver)\n")
} else {
  cat("RESULT: Significant difference — this contradicts your primer's earlier MAR finding, worth investigating\n")
}