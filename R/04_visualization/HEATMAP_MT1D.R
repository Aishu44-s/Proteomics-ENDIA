library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

q_threshold <- 0.10
number_to_plot <- 30

analysis_casecat <- "0"
target_term <- "corrected_visit_sample:maternal_T1D"

master_file <- paste0(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/",
  "NEW ANALYSIS/Helena_code/fixed_effects_table/",
  "fixed_effects_master_results_adj_int_rm.csv"
)

emm_file <- paste0(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/",
  "NEW ANALYSIS/Helena_code/emmeans_tables/",
  "emmeans_t1d_visit_casecat0_means_q10.csv"
)

out_dir <- paste0(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/",
  "NEW ANALYSIS/Helena_code/heatmap_plot"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# Select the top proteins from the original D3/D5 full model
# ---------------------------------------------------------------------------

fixed_effects_master <- read.csv(
  master_file,
  check.names = FALSE
)

selected_all <- fixed_effects_master %>%
  filter(
    Variable == target_term,
    ModelName == "form_adj_int_rm",
    is.finite(QValue),
    QValue < q_threshold
  ) %>%
  distinct(Protein, .keep_all = TRUE) %>%
  arrange(QValue)

if (nrow(selected_all) == 0) {
  stop(
    "No proteins pass q < 0.10 for the maternal-T1D-by-trimester ",
    "interaction in the full D3/D5 model."
  )
}

selected <- selected_all %>%
  slice_head(n = number_to_plot) %>%
  select(Protein, QValue)

message(
  "Selected ", nrow(selected),
  " proteins from the full D3/D5 model for the heatmap."
)

# ---------------------------------------------------------------------------
# Read casecat = 0 adjusted estimated means created by revised D6
# ---------------------------------------------------------------------------

emm_means <- read.csv(
  emm_file,
  check.names = FALSE
)

if (
  "casecat_condition" %in% names(emm_means) &&
  any(as.character(emm_means$casecat_condition) != analysis_casecat)
) {
  stop(
    "The D6 means file does not contain only the requested ",
    "casecat = ", analysis_casecat, " condition."
  )
}

heatmap_data <- emm_means %>%
  filter(
    Protein %in% selected$Protein,
    is.finite(emmean)
  ) %>%
  select(
    Protein,
    corrected_visit_sample,
    maternal_T1D,
    emmean
  ) %>%
  left_join(selected, by = "Protein")

if (nrow(heatmap_data) == 0) {
  stop(
    "No estimable casecat = 0 adjusted means were available ",
    "for the selected proteins."
  )
}

# ---------------------------------------------------------------------------
# Scale estimated means within each protein
# ---------------------------------------------------------------------------

heatmap_data <- heatmap_data %>%
  group_by(Protein) %>%
  mutate(
    MeanZ = if (n() > 1 && sd(emmean) > 0) {
      as.numeric(scale(emmean))
    } else {
      0
    }
  ) %>%
  ungroup() %>%
  mutate(
    Trimester = factor(
      corrected_visit_sample,
      levels = c("T01", "T02", "T03"),
      labels = c("Trimester 1", "Trimester 2", "Trimester 3")
    ),
    MaternalT1D = factor(
      as.character(maternal_T1D),
      levels = c("0", "1"),
      labels = c("No maternal T1D", "Maternal T1D")
    ),
    Group = factor(
      paste(Trimester, MaternalT1D, sep = "\n"),
      levels = c(
        "Trimester 1\nMaternal T1D",
        "Trimester 2\nMaternal T1D",
        "Trimester 3\nMaternal T1D",
        "Trimester 1\nNo maternal T1D",
        "Trimester 2\nNo maternal T1D",
        "Trimester 3\nNo maternal T1D"
      )
    ),
    GroupPosition = as.integer(Group)
  )

# Keep proteins ordered by q-value, with the smallest q-value at the top
protein_order <- selected %>%
  filter(Protein %in% heatmap_data$Protein) %>%
  arrange(desc(QValue)) %>%
  pull(Protein)

heatmap_data$Protein <- factor(
  heatmap_data$Protein,
  levels = protein_order
)

# ---------------------------------------------------------------------------
# Identify significant post-hoc maternal T1D vs no maternal T1D contrasts
# ---------------------------------------------------------------------------

posthoc_file <- paste0(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/",
  "NEW ANALYSIS/Helena_code/model_output_summary/",
  "maternal_T1D_by_trimester_q_lt_0.10_table.csv"
)

posthoc_results <- read.csv(posthoc_file, check.names = FALSE)

# Put the trimester-specific post-hoc q-values into long format
posthoc_stars <- bind_rows(
  posthoc_results %>%
    transmute(
      Protein,
      Trimester = "Trimester 1",
      PostHoc_QValue = Contrast_QValue_T1
    ),
  
  posthoc_results %>%
    transmute(
      Protein,
      Trimester = "Trimester 2",
      PostHoc_QValue = Contrast_QValue_T2
    ),
  
  posthoc_results %>%
    transmute(
      Protein,
      Trimester = "Trimester 3",
      PostHoc_QValue = Contrast_QValue_T3
    )
) %>%
  filter(
    Protein %in% protein_order,
    is.finite(PostHoc_QValue),
    PostHoc_QValue < q_threshold
  )

# One star per significant maternal-T1D versus no-maternal-T1D comparison.
# Place it only on the maternal-T1D tile to avoid repeating it twice.
star_cells <- heatmap_data %>%
  filter(MaternalT1D == "Maternal T1D") %>%
  left_join(
    posthoc_stars,
    by = c("Protein", "Trimester")
  ) %>%
  filter(!is.na(PostHoc_QValue))

# ---------------------------------------------------------------------------
# Create one combined heatmap:
# maternal T1D T1–T3, then no maternal T1D T1–T3
# ---------------------------------------------------------------------------

# Symmetric scale around zero: blue and red have equal visual meaning
z_max <- max(abs(heatmap_data$MeanZ), na.rm = TRUE)
z_limits <- c(-z_max, z_max)

p_combined_heatmap <- ggplot(
  heatmap_data,
  aes(x = Group, y = Protein, fill = MeanZ)
) +
  geom_tile(
    colour = "white",
    linewidth = 0.5,
    width = 0.95
  ) +
  
  geom_text(
    data = star_cells,
    aes(x = Group, y = Protein, label = "*"),
    inherit.aes = FALSE,
    colour = "black",
    fontface = "bold",
    size = 5
  ) +
  
  # Divider between maternal T1D and no maternal T1D columns
  geom_vline(
    xintercept = 3.5,
    colour = "white",
    linewidth = 2
  ) +
  
  scale_x_discrete(
    limits = levels(heatmap_data$Group),
    labels = c(
      "T1\nMaternal\nT1D",
      "T2\nMaternal\nT1D",
      "T3\nMaternal\nT1D",
      "T1\nNo maternal\nT1D",
      "T2\nNo maternal\nT1D",
      "T3\nNo maternal\nT1D"
    ),
    drop = FALSE
  ) +
  
  scale_y_discrete(drop = FALSE) +
  
  scale_fill_gradient2(
    low = "#2B6CB0",
    mid = "#FAFAF7",
    high = "#C53030",
    midpoint = 0,
    limits = z_limits,
    oob = scales::squish,
    name = "Within-protein\nz-score"
  ) +
  
  labs(
    title = "Adjusted protein trajectories across pregnancy",
    subtitle = paste0(
      "Top ", length(protein_order),
      " maternal-T1D-by-trimester interaction proteins",
      "\nEstimated within casecat = ", analysis_casecat
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "* FDR-adjusted post-hoc maternal T1D versus no maternal T1D ",
      "comparison at that trimester: q < 0.10."
    )
  ) +
  
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 11),
    axis.text.y = element_text(size = 13),
    plot.title = element_text(
      face = "bold",
      size = 18
    ),
    plot.subtitle = element_text(
      size = 12,
      lineheight = 1.1
    ) ,
    plot.caption = element_text(size = 10, hjust = 0),
    legend.position = "right"
  )

print(p_combined_heatmap)

# ---------------------------------------------------------------------------
# Save the combined heatmap
# ---------------------------------------------------------------------------

plot_height <- max(7, 0.32 * length(protein_order) + 3)

ggsave(
  filename = file.path(
    out_dir,
    paste0(
      "heatmap_combined_maternal_T1D_casecat",
      analysis_casecat,
      "_top",
      number_to_plot,
      ".png"
    )
  ),
  plot = p_combined_heatmap,
  width = 11,
  height = plot_height,
  units = "in",
  dpi = 300,
  bg = "white"
)