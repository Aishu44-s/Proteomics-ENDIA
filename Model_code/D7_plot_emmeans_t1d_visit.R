# ============================================================================
# D7_plot_emmeans_t1d_visit.R
# Run AFTER D6_bulk_emmeans_t1d_visit.R
# ============================================================================

library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

target_term <- "corrected_visit_sample:maternal_T1D"
q_threshold <- 0.05
number_of_top_proteins <- 12

plot_folder <- "plots_t1d_visit"

if (!dir.exists(plot_folder)) {
  dir.create(plot_folder)
}

# ---------------------------------------------------------------------------
# Read D5 and D6 output files
# ---------------------------------------------------------------------------

master_results <- read.csv(
  "fixed_effects_master_results_adj_int_rm.csv",
  check.names = FALSE
)

emm_means_long <- read.csv(
  "emmeans_t1d_visit_means_q05.csv",
  check.names = FALSE
)

# ---------------------------------------------------------------------------
# Identify significant proteins and retain their q-values
# ---------------------------------------------------------------------------

significant_proteins <- master_results %>%
  filter(
    Variable == target_term,
    QValue < q_threshold
  ) %>%
  select(Protein, QValue) %>%
  distinct() %>%
  arrange(QValue)

top_proteins <- significant_proteins %>%
  slice_head(n = number_of_top_proteins) %>%
  pull(Protein)

# Join q-values to the estimated marginal means
plot_data <- emm_means_long %>%
  inner_join(significant_proteins, by = "Protein") %>%
  mutate(
    maternal_T1D = factor(maternal_T1D),
    corrected_visit_sample = factor(corrected_visit_sample)
  )

# ============================================================================
# 1. INDIVIDUAL ADJUSTED TRAJECTORY PLOTS
# One plot per protein
# ============================================================================

for (protein_name in top_proteins) {
  
  one_protein <- plot_data %>%
    filter(Protein == protein_name)
  
  protein_q <- unique(one_protein$QValue)
  
  p <- ggplot(
    one_protein,
    aes(
      x = corrected_visit_sample,
      y = emmean,
      group = maternal_T1D,
      colour = maternal_T1D
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.5) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.10,
      linewidth = 0.6
    ) +
    labs(
      title = protein_name,
      subtitle = paste0(
        "Maternal T1D × trimester interaction q = ",
        signif(protein_q, 3)
      ),
      x = "Trimester",
      y = "Adjusted predicted protein abundance",
      colour = "Maternal T1D"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", protein_name)
  
  ggsave(
    filename = file.path(
      plot_folder,
      paste0("trajectory_", safe_name, ".png")
    ),
    plot = p,
    width = 7,
    height = 5,
    dpi = 300
  )
}

# ============================================================================
# 2. FACETED TRAJECTORY PLOT
# Top proteins together in one figure
# ============================================================================

top_plot_data <- plot_data %>%
  filter(Protein %in% top_proteins) %>%
  mutate(
    Protein = factor(Protein, levels = top_proteins)
  )

p_faceted <- ggplot(
  top_plot_data,
  aes(
    x = corrected_visit_sample,
    y = emmean,
    group = maternal_T1D,
    colour = maternal_T1D
  )
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.08,
    linewidth = 0.4
  ) +
  facet_wrap(~ Protein, scales = "free_y", ncol = 3) +
  labs(
    title = "Top proteins with a maternal-T1D-by-trimester interaction",
    x = "Trimester",
    y = "Adjusted predicted protein abundance",
    colour = "Maternal T1D"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    plot_folder,
    "top_proteins_t1d_by_trimester_faceted.png"
  ),
  plot = p_faceted,
  width = 12,
  height = 12,
  dpi = 300
)

# ============================================================================
# 3. HEATMAP
# Shows each protein's relative adjusted pattern across group x trimester
# ============================================================================

heatmap_data <- top_plot_data %>%
  mutate(
    Group_Trimester = paste(
      "T1D", maternal_T1D,
      "- Trimester", corrected_visit_sample
    )
  ) %>%
  group_by(Protein) %>%
  mutate(
    Z_score = as.numeric(scale(emmean))
  ) %>%
  ungroup() %>%
  mutate(
    Protein = factor(Protein, levels = rev(top_proteins))
  )

p_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = Group_Trimester,
    y = Protein,
    fill = Z_score
  )
) +
  geom_tile() +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Relative\nabundance\n(z-score)"
  ) +
  labs(
    title = "Adjusted protein patterns by maternal T1D status and trimester",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid = element_blank()
  )

ggsave(
  filename = file.path(
    plot_folder,
    "top_proteins_t1d_by_trimester_heatmap.png"
  ),
  plot = p_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)

# ============================================================================
# 4. HIT-COUNT BAR PLOT
# Number of q-significant proteins for every model term
# ============================================================================

hit_counts <- master_results %>%
  group_by(Variable) %>%
  summarise(
    SignificantProteins = sum(QValue < q_threshold, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(SignificantProteins)) %>%
  mutate(
    Variable = factor(Variable, levels = Variable)
  )

p_hits <- ggplot(
  hit_counts,
  aes(
    x = Variable,
    y = SignificantProteins
  )
) +
  geom_col(fill = "#2C7FB8") +
  coord_flip() +
  labs(
    title = paste0(
      "Proteins significant at q < ",
      q_threshold
    ),
    x = NULL,
    y = "Number of significant proteins"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(
    plot_folder,
    "significant_protein_hit_counts.png"
  ),
  plot = p_hits,
  width = 9,
  height = 6,
  dpi = 300
)