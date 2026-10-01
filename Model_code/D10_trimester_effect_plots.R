# ============================================================================
# D10_trimester_effect_plots.R
# Proteins that significantly change across pregnancy trimesters
# ============================================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(ggplot2)
library(tidyr)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

target_term <- "corrected_visit_sample"
q_threshold <- 0.10
number_of_top_proteins <- 12

model_name <- "form_adj_int_rm"
ci_level <- 0.90

plot_folder <- "plots_trimester_effect"

if (!dir.exists(plot_folder)) {
  dir.create(plot_folder)
}

# ---------------------------------------------------------------------------
# Read master table and select trimester-significant proteins
# ---------------------------------------------------------------------------

master_results <- read.csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Helena_code/fixed_effects_table/fixed_effects_master_results_adj_int_rm.csv",
  check.names = FALSE
)

selected_proteins <- master_results %>%
  filter(
    Variable == target_term,
    QValue < q_threshold
  ) %>%
  select(Protein, QValue) %>%
  distinct() %>%
  arrange(QValue)

if (nrow(selected_proteins) == 0) {
  stop("No proteins meet q < 0.10 for the trimester effect.")
}

top_proteins <- selected_proteins %>%
  slice_head(n = number_of_top_proteins) %>%
  pull(Protein)

message(
  "Number of proteins significant for trimester at q < 0.10: ",
  nrow(selected_proteins)
)

# ---------------------------------------------------------------------------
# Use the same adjusted model as D3/D5
# ---------------------------------------------------------------------------

if (!model_name %in% names(model)) {
  stop(paste("Could not find", model_name, "in model."))
}

model_formula <- as.formula(model[[model_name]])

# ---------------------------------------------------------------------------
# Helper for cleaner protein labels
# ---------------------------------------------------------------------------

prettify_protein_id <- function(x) {
  
  parts <- strsplit(x, "\\.")
  
  vapply(parts, function(p) {
    p <- p[p != "sp"]
    paste(p, collapse = " · ")
  }, character(1))
}

# ---------------------------------------------------------------------------
# Fit selected proteins and calculate overall trimester EMMs
# ---------------------------------------------------------------------------

emm_trimester_list <- list()
fit_log_list <- list()

for (protein in selected_proteins$Protein) {
  
  cat("Processing:", protein, "\n")
  
  protein_q <- selected_proteins %>%
    filter(Protein == protein) %>%
    pull(QValue)
  
  tmp <- mother
  tmp$btmean <- tmp[[protein]]
  tmp <- tmp[!is.na(tmp$btmean), ]
  
  fit <- tryCatch(
    lmer(
      formula = model_formula,
      data = tmp,
      weights = weightTEDDY,
      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
      )
    ),
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    
    fit_log_list[[protein]] <- data.frame(
      Protein = protein,
      Status = "Failed",
      Message = fit$message,
      stringsAsFactors = FALSE
    )
    
    next
  }
  
  fit_log_list[[protein]] <- data.frame(
    Protein = protein,
    Status = "Success",
    Message = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Overall adjusted mean at each trimester.
  # Proportional weights use the observed distribution of group variables.
  emm <- emmeans(
    fit,
    ~ corrected_visit_sample,
    weights = "proportional",
    lmer.df = "satterthwaite"
  )
  
  emm_df <- as.data.frame(
    summary(
      emm,
      infer = c(TRUE, TRUE),
      level = ci_level
    )
  )
  
  emm_df$Protein <- protein
  emm_df$TrimesterQValue <- protein_q
  
  emm_trimester_list[[protein]] <- emm_df
}

emm_trimester_long <- bind_rows(emm_trimester_list)
fit_log <- bind_rows(fit_log_list)

write.csv(
  emm_trimester_long,
  "emmeans_trimester_means_q10.csv",
  row.names = FALSE
)

write.csv(
  fit_log,
  "emmeans_trimester_fit_log_q10.csv",
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Prepare all significant trimester proteins for plotting
# ---------------------------------------------------------------------------

all_plot_data <- emm_trimester_long %>%
  mutate(
    corrected_visit_sample = factor(
      corrected_visit_sample,
      levels = c("T01", "T02", "T03"),
      labels = c("Trimester 1", "Trimester 2", "Trimester 3")
    ),
    
    ProteinLabel = prettify_protein_id(Protein)
  ) %>%
  filter(
    !is.na(emmean),
    !is.na(lower.CL),
    !is.na(upper.CL)
  )

# Keep the top 12 only for the main faceted figure
plot_data <- all_plot_data %>%
  filter(Protein %in% top_proteins) %>%
  mutate(
    Protein = factor(Protein, levels = top_proteins)
  )
# ---------------------------------------------------------------------------
# Shared plot theme
# ---------------------------------------------------------------------------

trimester_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(colour = "grey30"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#D9D9D9"),
    legend.position = "none"
  )

# ---------------------------------------------------------------------------
# Individual plots
# ---------------------------------------------------------------------------

for (protein_name in top_proteins) {
  
  one_protein <- plot_data %>%
    filter(Protein == protein_name)
  
  protein_q <- unique(one_protein$TrimesterQValue)
  protein_label <- unique(one_protein$ProteinLabel)
  
  p_individual <- ggplot(
    one_protein,
    aes(
      x = corrected_visit_sample,
      y = emmean,
      group = 1
    )
  ) +
    geom_errorbar(
      aes(ymin = lower.CL, ymax = upper.CL),
      width = 0.08,
      colour = "#0072B2",
      alpha = 0.55
    ) +
    geom_line(
      colour = "#0072B2",
      linewidth = 0.9
    ) +
    geom_point(
      colour = "#0072B2",
      size = 3
    ) +
    labs(
      title = protein_label,
      subtitle = paste0(
        "Overall trimester effect, q = ",
        format(signif(protein_q, 3), scientific = TRUE)
      ),
      x = NULL,
      y = "Adjusted predicted protein abundance"
    ) +
    trimester_theme
  
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", as.character(protein_name))
  
  ggsave(
    filename = file.path(
      plot_folder,
      paste0("trimester_trajectory_", safe_name, ".png")
    ),
    plot = p_individual,
    width = 6,
    height = 4.5,
    dpi = 300,
    bg = "white"
  )
}

# ---------------------------------------------------------------------------
# Faceted overview plot
# ---------------------------------------------------------------------------

protein_labels <- plot_data %>%
  distinct(Protein, ProteinLabel, TrimesterQValue) %>%
  mutate(
    FacetLabel = paste0(
      ProteinLabel,
      "\nq = ",
      format(TrimesterQValue, digits = 2, scientific = TRUE)
    )
  )

facet_labels <- setNames(
  protein_labels$FacetLabel,
  protein_labels$Protein
)

p_faceted <- ggplot(
  plot_data,
  aes(
    x = corrected_visit_sample,
    y = emmean,
    group = 1
  )
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.08,
    colour = "#0072B2",
    alpha = 0.50
  ) +
  geom_line(
    colour = "#0072B2",
    linewidth = 0.85
  ) +
  geom_point(
    colour = "#0072B2",
    size = 2.5
  ) +
  facet_wrap(
    ~ Protein,
    scales = "free_y",
    ncol = 3,
    labeller = labeller(Protein = facet_labels)
  ) +
  labs(
    title = "Top proteins with significant trimester effects",
    subtitle = "Adjusted estimated marginal means for estimable trimesters (T2 and T3), with 90% confidence intervals",
    x = NULL,
    y = "Adjusted predicted protein abundance"
  ) +
  trimester_theme

ggsave(
  filename = file.path(
    plot_folder,
    "top_proteins_trimester_effect_faceted.png"
  ),
  plot = p_faceted,
  width = 13,
  height = 14,
  dpi = 300,
  bg = "white"
)

# ============================================================================
# DESCRIPTIVE HEATMAP: ALL SIGNIFICANT TRIMESTER PROTEINS
# ============================================================================

# Use the objects already created earlier in your D10 script
sig_proteins <- selected_proteins$Protein

protein_q_table <- selected_proteins %>%
  select(Protein, QValue) %>%
  distinct()

heatmap_long <- mother %>%
  select(corrected_visit_sample, all_of(sig_proteins)) %>%
  filter(!is.na(corrected_visit_sample)) %>%
  mutate(
    corrected_visit_sample = factor(
      corrected_visit_sample,
      levels = c("T01", "T02", "T03"),
      labels = c("Trimester 1", "Trimester 2", "Trimester 3")
    )
  ) %>%
  pivot_longer(
    cols = all_of(sig_proteins),
    names_to = "Protein",
    values_to = "Abundance"
  ) %>%
  group_by(Protein, corrected_visit_sample) %>%
  summarise(
    MeanAbundance = mean(Abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(MeanAbundance)) %>%
  group_by(Protein) %>%
  filter(n() == 3) %>%  # retain proteins with T1, T2 and T3 values
  mutate(ZScore = as.numeric(scale(MeanAbundance))) %>%
  ungroup() %>%
  left_join(protein_q_table, by = "Protein")

# Cluster proteins with similar T1-to-T3 patterns
heatmap_matrix <- heatmap_long %>%
  select(Protein, corrected_visit_sample, ZScore) %>%
  pivot_wider(
    names_from = corrected_visit_sample,
    values_from = ZScore
  ) %>%
  arrange(Protein)

profile_matrix <- as.matrix(heatmap_matrix[, -1])
rownames(profile_matrix) <- heatmap_matrix$Protein

protein_order <- rownames(profile_matrix)[
  hclust(dist(profile_matrix), method = "complete")$order
]

# Save the row order, so each heatmap row can be identified
heatmap_order <- tibble(
  Protein = protein_order,
  HeatmapRow = seq_along(protein_order)
) %>%
  left_join(protein_q_table, by = "Protein")

write.csv(
  heatmap_order,
  file.path(plot_folder, "trimester_heatmap_protein_order_q10.csv"),
  row.names = FALSE
)

heatmap_plot_data <- heatmap_long %>%
  filter(Protein %in% protein_order) %>%
  mutate(
    Protein = factor(Protein, levels = rev(protein_order))
  )

p_heatmap_all <- ggplot(
  heatmap_plot_data,
  aes(
    x = corrected_visit_sample,
    y = Protein,
    fill = ZScore
  )
) +
  geom_tile(colour = "white", linewidth = 0.08) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "#F7F7F7",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    oob = scales::squish,
    name = "Relative protein\nabundance\n(z-score)"
  ) +
  labs(
    title = "Trimester patterns among significant maternal proteins",
    subtitle = paste0(
      "Observed mean abundance; ",
      nrow(heatmap_order),
      " proteins with adjusted trimester-effect q < ",
      q_threshold
    ),
    x = NULL,
    y = NULL,
    caption = paste(
      "Blue/red show abundance relative to each protein's own mean.",
      "Colours cannot be compared between different proteins."
    )
  ) +
  theme_minimal(base_size = 12) +
  scale_y_discrete(labels = prettify_protein_id) +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(
      size = 5,
      colour = "black",
      margin = margin(r = 4)
    ),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 13),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12),
    plot.caption = element_text(size = 10, colour = "grey35"),
    legend.position = "right",
    plot.margin = margin(10, 20, 10, 120, unit = "pt")
  )

ggsave(
  filename = file.path(
    plot_folder,
    "all_significant_trimester_proteins_descriptive_heatmap_q10_LABELLED.png"
  ),
  plot = p_heatmap_all,
  width = 16,
  height = 48,
  dpi = 300,
  bg = "white",
  limitsize = FALSE
)
