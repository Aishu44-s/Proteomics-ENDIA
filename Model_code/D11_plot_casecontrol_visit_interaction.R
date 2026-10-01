# ============================================================================
# D11_plot_casecontrol_visit_interaction.R
# Do case and control infants show different maternal protein trajectories?
# ============================================================================

library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(ggplot2)

# ---------------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------------

target_term <- "casecat:corrected_visit_sample"
q_threshold <- 0.10
number_of_top_proteins <- 12

# Use the same complete adjusted model used for your master table
model_name <- "form_adj_int_rm"
ci_level <- 0.90

plot_folder <- "plots_casecontrol_trajectory"

if (!dir.exists(plot_folder)) {
  dir.create(plot_folder)
}

# ---------------------------------------------------------------------------
# Read master table and identify significant case/control-by-trimester proteins
# ---------------------------------------------------------------------------

master_results <- read.csv(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Helena_code/fixed_effects_table/fixed_effects_master_results_adj_int_rm.csv",
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

message(
  "Number of proteins significant for the case/control-by-trimester interaction ",
  "at q < ", q_threshold, ": ", nrow(selected_proteins)
)

if (nrow(selected_proteins) == 0) {
  stop(
    "No proteins meet q < ", q_threshold,
    " for the case/control-by-trimester interaction. ",
    "Do not create case/control trajectory plots as a primary result."
  )
}

top_proteins <- selected_proteins %>%
  slice_head(n = number_of_top_proteins) %>%
  pull(Protein)

# ---------------------------------------------------------------------------
# Retrieve the model formula used in D3
# ---------------------------------------------------------------------------

if (!model_name %in% names(model)) {
  stop(paste("Could not find", model_name, "in model."))
}

model_formula <- as.formula(model[[model_name]])

# ---------------------------------------------------------------------------
# Fit significant proteins and obtain adjusted means
# ---------------------------------------------------------------------------

emm_means_list <- list()
emm_contrasts_list <- list()
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
      Message = fit$message
    )
    
    next
  }
  
  fit_log_list[[protein]] <- data.frame(
    Protein = protein,
    Status = "Success",
    Message = NA_character_
  )
  
  # Adjusted mean at each trimester, separately for case/control infants
  emm <- emmeans(
    fit,
    ~ corrected_visit_sample | casecat,
    lmer.df = "satterthwaite"
  )
  
  means_df <- as.data.frame(
    summary(
      emm,
      infer = c(TRUE, TRUE),
      level = ci_level
    )
  )
  
  means_df$Protein <- protein
  means_df$CaseVisitInteractionQValue <- protein_q
  
  # Within each infant group: T1 vs T2, T1 vs T3, T2 vs T3
  contrasts_df <- as.data.frame(
    summary(
      pairs(emm, adjust = "tukey"),
      infer = c(TRUE, TRUE),
      level = ci_level
    )
  )
  
  contrasts_df$Protein <- protein
  contrasts_df$CaseVisitInteractionQValue <- protein_q
  
  emm_means_list[[protein]] <- means_df
  emm_contrasts_list[[protein]] <- contrasts_df
}

emm_means <- bind_rows(emm_means_list)
emm_contrasts <- bind_rows(emm_contrasts_list)
fit_log <- bind_rows(fit_log_list)

write.csv(
  emm_means,
  file.path(plot_folder, "emmeans_casecontrol_visit_means_q10.csv"),
  row.names = FALSE
)

write.csv(
  emm_contrasts,
  file.path(plot_folder, "emmeans_casecontrol_visit_pairwise_contrasts_q10.csv"),
  row.names = FALSE
)

write.csv(
  fit_log,
  file.path(plot_folder, "emmeans_casecontrol_visit_fit_log_q10.csv"),
  row.names = FALSE
)

# ---------------------------------------------------------------------------
# Prepare data for plotting
# ---------------------------------------------------------------------------

plot_data <- emm_means %>%
  filter(Protein %in% top_proteins) %>%
  mutate(
    corrected_visit_sample = factor(
      corrected_visit_sample,
      levels = c("T01", "T02", "T03"),
      labels = c("Trimester 1", "Trimester 2", "Trimester 3")
    ),
    casecat = factor(
      as.character(casecat),
      levels = c("0", "1"),
      labels = c("Control infant", "Case infant")
    )
  ) %>%
  filter(
    !is.na(emmean),
    !is.na(lower.CL),
    !is.na(upper.CL)
  )

if (nrow(plot_data) == 0) {
  stop("No estimable adjusted means were available for plotting.")
}

protein_labels <- plot_data %>%
  distinct(Protein, CaseVisitInteractionQValue) %>%
  mutate(
    ProteinLabel = gsub("^sp\\.", "", Protein),
    FacetLabel = paste0(
      ProteinLabel,
      "\nq = ",
      format(CaseVisitInteractionQValue, digits = 2, scientific = TRUE)
    )
  )

facet_labels <- setNames(
  protein_labels$FacetLabel,
  protein_labels$Protein
)

# ---------------------------------------------------------------------------
# Faceted case/control trajectory plot
# ---------------------------------------------------------------------------

case_colours <- c(
  "Control infant" = "#0072B2",
  "Case infant" = "#D55E00"
)

case_shapes <- c(
  "Control infant" = 16,
  "Case infant" = 17
)

case_linetypes <- c(
  "Control infant" = "solid",
  "Case infant" = "dashed"
)

p_casecontrol <- ggplot(
  plot_data,
  aes(
    x = corrected_visit_sample,
    y = emmean,
    colour = casecat,
    shape = casecat,
    linetype = casecat,
    group = casecat
  )
) +
  geom_errorbar(
    aes(ymin = lower.CL, ymax = upper.CL),
    width = 0.08,
    alpha = 0.45
  ) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.8) +
  facet_wrap(
    ~ Protein,
    scales = "free_y",
    ncol = 3,
    labeller = labeller(Protein = facet_labels)
  ) +
  scale_colour_manual(values = case_colours) +
  scale_shape_manual(values = case_shapes) +
  scale_linetype_manual(values = case_linetypes) +
  labs(
    title = "Maternal protein trajectories by infant case/control status",
    subtitle = "Adjusted estimated marginal means with 90% confidence intervals",
    x = "Trimester",
    y = "Adjusted predicted protein abundance",
    colour = "Infant group",
    shape = "Infant group",
    linetype = "Infant group"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(colour = "grey30"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "#F2F2F2", colour = NA),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#D9D9D9"),
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    plot_folder,
    "top_proteins_casecontrol_by_trimester_faceted_q10.png"
  ),
  plot = p_casecontrol,
  width = 14,
  height = 15,
  dpi = 300,
  bg = "white"
)