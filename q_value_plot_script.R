library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

# ---- 0. Settings ----------------------------------------------------------

selected_model <- "form_adj_int_rm"   # change this if you re-run for a different tier
q_thresh_strict <- 0.05
q_thresh_loose  <- 0.10

qval_files <- c(
  cc         = "qvalue_cc_cal.csv",
  visit      = "qvalue_visit_cal.csv",
  ccvisit    = "qvalue_ccvisit_cal.csv",
  t1dcc      = "qvalue_t1dcc_cal.csv",
  t1dvisit   = "qvalue_t1dvisit_cal.csv",
  t1dccvisit = "qvalue_t1dccvisit_cal.csv"
)

dir.create("qvalue_plots", showWarnings = FALSE)

# ---- 1. Read the selected model's row out of each q-value matrix ----------

read_qrow <- function(file, term) {
  df <- read_csv(file, show_col_types = FALSE)
  names(df)[1] <- "model"
  row <- df %>% filter(model == selected_model)
  if (nrow(row) == 0) {
    stop(sprintf("Model '%s' not found in %s -- check rownames.", selected_model, file))
  }
  protein <- names(row)[-1]
  qvalue  <- as.numeric(row[1, -1])
  tibble(protein = protein, term = term, qvalue = qvalue)
}

qlong <- bind_rows(
  lapply(names(qval_files), function(term) read_qrow(qval_files[[term]], term))
)

qlong <- qlong %>%
  mutate(
    sig_05 = qvalue < q_thresh_strict,
    sig_10 = qvalue < q_thresh_loose
  )

# Long format -- one row per protein x term. Easiest to filter/sort.
write_csv(qlong, sprintf("qvalue_long_%s.csv", selected_model))

# Wide format -- one row per protein, one column per term. Easier to eyeball
# a specific protein across all six terms at once.
qwide <- qlong %>%
  select(protein, term, qvalue) %>%
  pivot_wider(names_from = term, values_from = qvalue)
write_csv(qwide, sprintf("qvalue_wide_%s.csv", selected_model))

# ---- 2. Hit counts per term ------------------------------------------------

summary_counts <- qlong %>%
  group_by(term) %>%
  summarise(
    n_sig_05 = sum(sig_05, na.rm = TRUE),
    n_sig_10 = sum(sig_10, na.rm = TRUE),
    .groups = "drop"
  )
write_csv(summary_counts, "qvalue_hit_counts.csv")
print(summary_counts)

p_counts <- summary_counts %>%
  pivot_longer(cols = c(n_sig_05, n_sig_10), names_to = "threshold", values_to = "n") %>%
  mutate(threshold = recode(threshold, n_sig_05 = "q < 0.05", n_sig_10 = "q < 0.10")) %>%
  ggplot(aes(x = term, y = n, fill = threshold)) +
  geom_col(position = "dodge") +
  labs(title = sprintf("Significant proteins per term (%s)", selected_model),
       x = "Term", y = "Number of proteins", fill = "Threshold") +
  theme_minimal()
ggsave("qvalue_plots/hit_counts_barplot.jpeg", p_counts, width = 7, height = 5, dpi = 200)

# ---- 3. Ranked q-value plot per term (a quick "how many, how strong" view) -

qlong_ranked <- qlong %>%
  group_by(term) %>%
  arrange(qvalue) %>%
  mutate(rank = row_number()) %>%
  ungroup()

p_ranked <- ggplot(qlong_ranked, aes(x = rank, y = -log10(qvalue), color = sig_10)) +
  geom_point(size = 0.8) +
  geom_hline(yintercept = -log10(q_thresh_loose), linetype = "dashed") +
  geom_hline(yintercept = -log10(q_thresh_strict), linetype = "dotted") +
  facet_wrap(~term, scales = "free_x") +
  scale_color_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey60"),
                     labels = c(`TRUE` = "q<0.10", `FALSE` = "n.s."), name = NULL) +
  labs(title = sprintf("Ranked -log10(q) per term (%s)", selected_model),
       x = "Protein rank (sorted by q-value, within term)", y = "-log10(q-value)") +
  theme_minimal()
ggsave("qvalue_plots/ranked_qvalue_plot.jpeg", p_ranked, width = 10, height = 7, dpi = 200)

# ---- 4. Overlap between the two "live" terms: visit vs t1dvisit -----------
# (cc, ccvisit, t1dcc, t1dccvisit all had 0 hits for form_adj_int_rm, so the
#  only overlap worth checking right now is visit vs t1dvisit.)

visit_hits    <- qlong %>% filter(term == "visit",    sig_10) %>% pull(protein)
t1dvisit_hits <- qlong %>% filter(term == "t1dvisit", sig_10) %>% pull(protein)

overlap       <- intersect(visit_hits, t1dvisit_hits)
only_visit    <- setdiff(visit_hits, t1dvisit_hits)
only_t1dvisit <- setdiff(t1dvisit_hits, visit_hits)

cat("\nVisit-only hits:   ", length(only_visit), "\n")
cat("T1Dvisit-only hits:", length(only_t1dvisit), "\n")
cat("Overlap:           ", length(overlap), "\n")

overlap_summary <- tibble(
  category = c("visit only", "t1dvisit only", "overlap"),
  n = c(length(only_visit), length(only_t1dvisit), length(overlap))
)
write_csv(overlap_summary, "qvalue_plots/visit_t1dvisit_overlap_counts.csv")

p_overlap <- ggplot(overlap_summary, aes(x = category, y = n, fill = category)) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n), vjust = -0.3) +
  labs(title = "Overlap: visit hits vs t1dvisit hits (q<0.10)", x = NULL, y = "Number of proteins") +
  theme_minimal()
ggsave("qvalue_plots/visit_t1dvisit_overlap_barplot.jpeg", p_overlap, width = 6, height = 5, dpi = 200)

overlap_df <- tibble(
  protein = c(overlap, only_visit, only_t1dvisit),
  category = c(rep("overlap", length(overlap)),
               rep("visit_only", length(only_visit)),
               rep("t1dvisit_only", length(only_t1dvisit)))
)
write_csv(overlap_df, "qvalue_plots/visit_t1dvisit_overlap_proteins.csv")

cat("\nDone. Outputs:\n",
    "- qvalue_long_", selected_model, ".csv        (protein x term, tidy long format)\n",
    "- qvalue_wide_", selected_model, ".csv        (protein x term, wide format)\n",
    "- qvalue_hit_counts.csv                        (n significant per term, two thresholds)\n",
    "- qvalue_plots/hit_counts_barplot.jpeg\n",
    "- qvalue_plots/ranked_qvalue_plot.jpeg\n",
    "- qvalue_plots/visit_t1dvisit_overlap_barplot.jpeg\n",
    "- qvalue_plots/visit_t1dvisit_overlap_proteins.csv  (which proteins are in which bucket)\n",
    sep = "")