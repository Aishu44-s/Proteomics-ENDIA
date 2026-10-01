# ============================================================================
# D9_casecontrol_summary_plot.R
# Question 1:
# Are maternal proteins associated with infant case/control status?
# ============================================================================

library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

q_threshold <- 0.10
number_to_show <- 20

plot_folder <- "plots_casecontrol"

if (!dir.exists(plot_folder)) {
  dir.create(plot_folder)
}

# ---------------------------------------------------------------------------
# Read master results table
# ---------------------------------------------------------------------------

master_results <- read.csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Helena_code/fixed_effects_table/fixed_effects_master_results_adj_int_rm.csv",
  check.names = FALSE
)

# ---------------------------------------------------------------------------
# Extract the case/control results
# ---------------------------------------------------------------------------

casecat_results <- master_results %>%
  filter(Variable == "casecat") %>%
  mutate(
    Significant_q10 = QValue < q_threshold
  ) %>%
  arrange(QValue)

number_tested <- nrow(casecat_results)

number_significant <- casecat_results %>%
  filter(Significant_q10) %>%
  nrow()

message("Proteins tested for case/control: ", number_tested)
message("Proteins significant at q < 0.10: ", number_significant)

# Save all case/control results
write.csv(
  casecat_results,
  "casecontrol_all_protein_results.csv",
  row.names = FALSE
)

# Save q < 0.10 case/control hits only
casecat_hits_q10 <- casecat_results %>%
  filter(Significant_q10)

write.csv(
  casecat_hits_q10,
  "casecontrol_significant_proteins_q10.csv",
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Plot 1: case/control hit-count summary
# ---------------------------------------------------------------------------

hit_summary <- data.frame(
  Result = "Case/control",
  SignificantProteins = number_significant
)

p_hit_count <- ggplot(
  hit_summary,
  aes(x = Result, y = SignificantProteins)
) +
  geom_col(
    fill = "#0072B2",
    width = 0.55
  ) +
  geom_text(
    aes(label = SignificantProteins),
    vjust = -0.4,
    size = 5
  ) +
  labs(
    title = "Maternal proteins associated with infant case/control status",
    subtitle = paste0(
      "FDR threshold: q < ", q_threshold,
      "; proteins tested: ", number_tested
    ),
    x = NULL,
    y = "Number of significant proteins"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = file.path(
    plot_folder,
    "casecontrol_significant_protein_count_q10.png"
  ),
  plot = p_hit_count,
  width = 7,
  height = 5,
  dpi = 300,
  bg = "white"
)

# ---------------------------------------------------------------------------
# Plot 2: top 20 lowest-q case/control proteins
# This is descriptive. It shows the closest proteins to significance.
# ---------------------------------------------------------------------------

top_casecat <- casecat_results %>%
  slice_head(n = number_to_show) %>%
  mutate(
    Protein = factor(Protein, levels = rev(Protein)),
    MinusLog10Q = -log10(QValue),
    ResultLabel = ifelse(
      Significant_q10,
      "q < 0.10",
      "Not q-significant"
    )
  )

p_top_qvalues <- ggplot(
  top_casecat,
  aes(
    x = MinusLog10Q,
    y = Protein,
    colour = ResultLabel
  )
) +
  geom_point(size = 3) +
  geom_vline(
    xintercept = -log10(q_threshold),
    linetype = "dashed",
    colour = "#D55E00"
  ) +
  scale_colour_manual(
    values = c(
      "q < 0.10" = "#D55E00",
      "Not q-significant" = "#666666"
    )
  ) +
  labs(
    title = "Top case/control-associated proteins",
    subtitle = paste0(
      "Dashed line = q < ", q_threshold,
      " threshold"
    ),
    x = expression(-log[10](q-value)),
    y = NULL,
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    plot_folder,
    "casecontrol_top20_qvalues.png"
  ),
  plot = p_top_qvalues,
  width = 9,
  height = 7,
  dpi = 300,
  bg = "white"
)