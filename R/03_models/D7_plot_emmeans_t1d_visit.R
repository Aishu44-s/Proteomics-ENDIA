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
q_threshold <- 0.10
number_of_top_proteins <- 12

plot_folder <- "plots_t1d_visit_casecat0"

if (!dir.exists(plot_folder)) {
  dir.create(plot_folder)
}

# Colour-blind-friendly scheme
t1d_colours <- c(
  "No maternal T1D" = "#91bfdb",
  "Maternal T1D" = "#fc8d59"
)

t1d_shapes <- c(
  "No maternal T1D" = 16,
  "Maternal T1D" = 17
)

t1d_linetypes <- c(
  "No maternal T1D" = "solid",
  "Maternal T1D" = "dashed"
)

# ---------------------------------------------------------------------------
# Read results created by D5 and D6
# ---------------------------------------------------------------------------

master_results <- read.csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Helena_code/fixed_effects_table/fixed_effects_master_results_adj_int_rm.csv",
  check.names = FALSE
)

emm_means_long <- read.csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Helena_code/emmeans_tables/emmeans_t1d_visit_casecat0_means_q10.csv",
  check.names = FALSE
)

# ---------------------------------------------------------------------------
# Select significant proteins for the chosen interaction
# ---------------------------------------------------------------------------

significant_proteins <- master_results %>%
  filter(
    Variable == target_term,
    QValue < q_threshold
  ) %>%
  select(Protein, QValue) %>%
  distinct() %>%
  arrange(QValue)

if (nrow(significant_proteins) == 0) {
  stop("No proteins meet the selected q-value threshold.")
}

top_proteins <- significant_proteins %>%
  slice_head(n = number_of_top_proteins) %>%
  pull(Protein)

# ---------------------------------------------------------------------------
# Protein-label helper
# Example: sp.A0A075B6I0.LV861 becomes A0A075B6I0 · LV861
# ---------------------------------------------------------------------------

prettify_protein_id <- function(x) {
  
  parts <- strsplit(x, "\\.")
  
  vapply(parts, function(p) {
    p <- p[p != "sp"]
    paste(p, collapse = " · ")
  }, character(1))
}

# ---------------------------------------------------------------------------
# Prepare plotting data
# ---------------------------------------------------------------------------

plot_data <- emm_means_long %>%
  inner_join(significant_proteins, by = "Protein") %>%
  mutate(
    maternal_T1D = recode(
      as.character(maternal_T1D),
      "0" = "No maternal T1D",
      "1" = "Maternal T1D"
    ),
    
    maternal_T1D = factor(
      maternal_T1D,
      levels = c("No maternal T1D", "Maternal T1D")
    ),
    
    corrected_visit_sample = factor(
      corrected_visit_sample,
      levels = c("T01", "T02", "T03"),
      labels = c("Trimester 1", "Trimester 2", "Trimester 3")
    ),
    
    ProteinLabel = prettify_protein_id(Protein)
  )

# Check whether any estimated marginal means are unavailable
missing_emmeans <- plot_data %>%
  filter(is.na(emmean)) %>%
  distinct(maternal_T1D, corrected_visit_sample)

if (nrow(missing_emmeans) > 0) {
  
  message(
    "Some adjusted means are not estimable and will not be plotted:\n",
    paste(
      paste(
        missing_emmeans$maternal_T1D,
        "-",
        missing_emmeans$corrected_visit_sample
      ),
      collapse = "\n"
    )
  )
}

# Keep only estimable means for plotting
plot_data_usable <- plot_data %>%
  filter(
    !is.na(emmean),
    !is.na(lower.CL),
    !is.na(upper.CL)
  )

# ---------------------------------------------------------------------------
# Shared plot theme
# ---------------------------------------------------------------------------

theme_t1d_trajectory <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 15,
      margin = margin(b = 4)
    ),
    
    plot.subtitle = element_text(
      colour = "grey30",
      margin = margin(b = 12)
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 10,
      lineheight = 0.95
    ),
    
    strip.background = element_rect(
      fill = "#F2F2F2",
      colour = NA
    ),
    
    panel.grid.minor = element_blank(),
    
    panel.grid.major.x = element_blank(),
    
    panel.grid.major.y = element_line(
      colour = "#D9D9D9",
      linewidth = 0.3
    ),
    
    panel.spacing = unit(1.1, "lines"),
    
    legend.position = "top",
    
    legend.text = element_text(size = 11),
    
    axis.text.x = element_text(
      face = "bold",
      colour = "grey20"
    ),
    
    axis.title.y = element_text(
      margin = margin(r = 10)
    )
  )

# ============================================================================
# 1. Individual trajectory plot for each top protein
# ============================================================================

for (protein_name in top_proteins) {
  
  one_protein <- plot_data_usable %>%
    filter(Protein == protein_name)
  
  protein_q <- significant_proteins %>%
    filter(Protein == protein_name) %>%
    pull(QValue)
  
  display_name <- unique(one_protein$ProteinLabel)
  
  p_individual <- ggplot(
    one_protein,
    aes(
      x = corrected_visit_sample,
      y = emmean,
      group = maternal_T1D,
      colour = maternal_T1D,
      shape = maternal_T1D,
      linetype = maternal_T1D
    )
  ) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.08,
      linewidth = 0.5,
      alpha = 0.50
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 3) +
    scale_colour_manual(values = t1d_colours) +
    scale_shape_manual(values = t1d_shapes) +
    scale_linetype_manual(values = t1d_linetypes) +
    labs(
      title = display_name,
      subtitle = paste0(
        "Maternal T1D × trimester interaction, q = ",
        format(signif(protein_q, 3), scientific = TRUE)
      ),
      x = NULL,
      y = "Adjusted predicted protein abundance",
      colour = NULL,
      shape = NULL,
      linetype = NULL
    ) +
    theme_t1d_trajectory
  
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", protein_name)
  
  ggsave(
    filename = file.path(
      plot_folder,
      paste0("trajectory_", safe_name, ".png")
    ),
    plot = p_individual,
    width = 6,
    height = 4.5,
    dpi = 300,
    bg = "white"
  )
}

# ============================================================================
# 2. Faceted overview plot
# ============================================================================

top_plot_data <- plot_data_usable %>%
  filter(Protein %in% top_proteins) %>%
  mutate(
    Protein = factor(Protein, levels = top_proteins)
  )

protein_q_labels <- significant_proteins %>%
  filter(Protein %in% top_proteins) %>%
  mutate(
    ProteinLabel = prettify_protein_id(Protein),
    FacetLabel = paste0(
      ProteinLabel,
      "\nq = ",
      format(QValue, digits = 2, scientific = TRUE)
    )
  )

facet_labels <- setNames(
  protein_q_labels$FacetLabel,
  protein_q_labels$Protein
)

p_faceted_pretty <- ggplot(
  top_plot_data,
  aes(
    x = corrected_visit_sample,
    y = emmean,
    group = maternal_T1D,
    colour = maternal_T1D,
    shape = maternal_T1D,
    linetype = maternal_T1D
  )
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.08,
    linewidth = 0.45,
    alpha = 0.45
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  facet_wrap(
    ~ Protein,
    scales = "free_y",
    ncol = 3,
    labeller = labeller(Protein = facet_labels)
  ) +
  scale_colour_manual(values = t1d_colours) +
  scale_shape_manual(values = t1d_shapes) +
  scale_linetype_manual(values = t1d_linetypes) +
  labs(
    title = "Proteins with significant maternal-T1D-by-trimester interactions",
    subtitle = paste(
      "Adjusted estimated marginal means with 90% confidence intervals.",
      "Estimated within casecat = 0."
    ),
    x = NULL,
    y = "Adjusted predicted protein abundance",
    colour = NULL,
    shape = NULL,
    linetype = NULL
  ) +
  theme_t1d_trajectory

ggsave(
  filename = file.path(
    plot_folder,
    "top_proteins_t1d_by_trimester_casecat0.png"
  ),
  plot = p_faceted_pretty,
  width = 13,
  height = 14,
  dpi = 300,
  bg = "white"
)