library(dplyr)
library(readr)
library(mixOmics)
library(vegan)
library(ggplot2)

# ---- 1. Load data -------------------------------------------------------------

mother <- read_csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/protein_normalised_VSN.csv")

protein_col_pattern <- "^(sp|tr)\\."
resp <- names(mother)[grepl(protein_col_pattern, names(mother))]
cat("Number of proteins:", length(resp), "\n")

mother$corrected_visit_sample <- factor(mother$corrected_visit_sample)
mother$casecat <- factor(mother$case, levels = c(0,1), labels = c("Control","Case"))
# NOTE: keep this label distinct from maternal_T1D's labels below -- these
# are two independent variables (the label-mixup bug already fixed
# elsewhere in this project).
mother$T1Dcat  <- factor(mother$maternal_T1D, levels = c(0,1), labels = c("No_T1D","T1D"))

dir.create("multivariate_output", showWarnings = FALSE)

# ---- 2. Reusable function: PLS-DA + PERMANOVA + Christmas tree + top 20 -----

run_plsda_permanova <- function(data, group_col, group_label, trimester_label,
                                min_group_n = 5) {
  
  cat("\n=====================================================\n")
  cat("Trimester:", trimester_label, "| Grouping variable:", group_label, "\n")
  cat("=====================================================\n")
  
  Y <- data[[group_col]]
  X <- as.matrix(data[, resp])
  
  tab <- table(Y)
  cat("Group sizes:\n"); print(tab)
  
  if (any(tab < min_group_n)) {
    cat(">> SKIPPED: at least one group has fewer than", min_group_n,
        "samples -- PLS-DA/cross-validation would not be meaningful here.\n")
    return(NULL)
  }
  
  # ---- PLS-DA ----
  # ncomp upper bound: keep modest given small per-trimester sample sizes
  max_ncomp <- min(5, nrow(X) - 1)
  
  # Cross-validation folds: capped by the smallest group size so no fold
  # ends up empty for the minority class. If the minority group is very
  # small, fall back to a fixed ncomp=2 rather than attempting CV.
  min_grp_n <- min(tab)
  if (min_grp_n >= 4) {
    folds <- min(5, min_grp_n)
    plsda_init <- plsda(X, Y, ncomp = max_ncomp)
    set.seed(123)
    perf_res <- perf(plsda_init, validation = "Mfold", folds = folds,
                     nrepeat = 50, progressBar = FALSE)
    # choose ncomp minimising balanced error rate (BER) on the max-distance metric
    ncomp_chosen <- which.min(perf_res$error.rate$BER[, "max.dist"])
    cat("Cross-validated ncomp chosen:", ncomp_chosen, "\n")
  } else {
    ncomp_chosen <- 2
    cat("Group too small for reliable CV -- using fixed ncomp = 2.\n")
  }
  
  plsda_res <- plsda(X, Y, ncomp = ncomp_chosen)
  
  # ---- plotIndiv (sample scores/ordination) ----
  jpeg(sprintf("multivariate_output/plsda_indiv_%s_%s.jpeg", trimester_label, group_label),
       width = 1600, height = 1400, res = 200)
  plotIndiv(plsda_res, comp = c(1, min(2, ncomp_chosen)), group = Y,
            ind.names = FALSE,
            ellipse = TRUE, legend = TRUE,
            title = sprintf("PLS-DA: %s (%s)", group_label, trimester_label))
  dev.off()
  
  # ---- Christmas tree plot (plotLoadings, comp 1, top 20) ----
  jpeg(sprintf("multivariate_output/christmastree_%s_%s.jpeg", trimester_label, group_label),
       width = 1800, height = 2000, res = 200)
  plotLoadings(plsda_res, comp = 1, contrib = "max", method = "median",
               ndisplay = 20, size.name = 0.6,
               title = sprintf("Top 20 contributing proteins: %s (%s)",
                               group_label, trimester_label))
  dev.off()
  
  # ---- Top 20 contributors as a table (not just the plot) ----
  loadings_comp1 <- selectVar(plsda_res, comp = 1)$value
  loadings_comp1$protein <- rownames(loadings_comp1)
  top20 <- loadings_comp1 %>%
    arrange(desc(abs(value.var))) %>%
    slice_head(n = 20) %>%
    dplyr::select(protein, loading = value.var)
  
  write_csv(top20,
            sprintf("multivariate_output/top20_%s_%s.csv", trimester_label, group_label))
  
  # ---- PERMANOVA ----
  # Euclidean distance on the VSN-normalised protein matrix. (This is
  # continuous targeted-proteomics intensity data, not compositional
  # relative-abundance data, so plain Euclidean is appropriate here rather
  # than a compositional/Aitchison distance.)
  dist_mat <- dist(X, method = "euclidean")
  set.seed(123)
  permanova_res <- adonis2(dist_mat ~ Y, method = "euclidean", permutations = 999)
  print(permanova_res)
  
  # ---- PCoA (Principal Coordinates Analysis) ----
  # An UNSUPERVISED ordination on the same distance matrix used for PERMANOVA.
  # PLS-DA above is supervised -- it is fit knowing the group labels, and
  # will lean towards separating groups even when the underlying separation
  # is weak. PCoA does not use the group labels to build the axes at all; it
  # just lays out samples by overall similarity, and group labels are only
  # added afterwards as colour. If groups separate on a PCoA plot, that is a
  # much more convincing, label-independent confirmation of what PLS-DA and
  # PERMANOVA are suggesting. If PCoA shows total overlap while PLS-DA shows
  # separation, treat the PLS-DA separation with more caution.
  n_axes <- min(3, nrow(X) - 1)
  pcoa_res <- cmdscale(dist_mat, k = n_axes, eig = TRUE)
  eig_pos <- pcoa_res$eig[pcoa_res$eig > 0]
  var_explained <- round(100 * pcoa_res$eig[1:2] / sum(eig_pos), 1)
  
  pcoa_df <- data.frame(
    PCo1  = pcoa_res$points[, 1],
    PCo2  = pcoa_res$points[, 2],
    Group = Y
  )
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PCo1, y = PCo2, color = Group)) +
    geom_point(size = 2) +
    stat_ellipse(level = 0.68) +
    labs(title = sprintf("PCoA: %s (%s)", group_label, trimester_label),
         x = sprintf("PCo1 (%.1f%%)", var_explained[1]),
         y = sprintf("PCo2 (%.1f%%)", var_explained[2])) +
    theme_minimal()
  ggsave(sprintf("multivariate_output/pcoa_%s_%s.jpeg", trimester_label, group_label),
         p_pcoa, width = 7, height = 6, dpi = 200)
  
  # ---- PERMDISP (homogeneity of multivariate dispersion) ----
  # PERMANOVA (adonis2 above) tests whether group CENTROIDS differ. It
  # assumes the groups have similar within-group SPREAD around those
  # centroids. If that assumption is violated, a "significant" PERMANOVA
  # result can actually be driven by one group just being more variable
  # (more spread out) than the other, rather than by a genuine difference in
  # average protein profile. PERMDISP (vegan::betadisper + permutest) is the
  # multivariate analogue of Levene's test: it tests whether the average
  # distance of samples to their own group centroid differs between groups.
  #   - PERMDISP p > 0.05: dispersions are similar -- your PERMANOVA result
  #     can be interpreted as a genuine location (centroid) difference.
  #   - PERMDISP p < 0.05: dispersions differ -- a significant PERMANOVA
  #     result here is ambiguous (could be location, could be spread, could
  #     be both), and this should be flagged alongside the PERMANOVA result
  #     rather than silently ignored.
  disp_res <- betadisper(dist_mat, Y)
  set.seed(123)
  permdisp_res <- permutest(disp_res, permutations = 999)
  print(permdisp_res)
  
  jpeg(sprintf("multivariate_output/permdisp_boxplot_%s_%s.jpeg", trimester_label, group_label),
       width = 1200, height = 1000, res = 200)
  boxplot(disp_res, main = sprintf("PERMDISP: %s (%s)", group_label, trimester_label),
          ylab = "Distance to group centroid")
  dev.off()
  
  # ---- Return a one-row summary ----
  data.frame(
    trimester    = trimester_label,
    group_var    = group_label,
    n_total      = nrow(X),
    group_sizes  = paste(names(tab), tab, sep = "=", collapse = "; "),
    ncomp        = ncomp_chosen,
    permanova_R2 = permanova_res$R2[1],
    permanova_p  = permanova_res$`Pr(>F)`[1],
    pcoa_var1    = var_explained[1],
    pcoa_var2    = var_explained[2],
    permdisp_p   = permdisp_res$tab$`Pr(>F)`[1]
  )
}

# ---- 3. Run across all trimester x grouping-variable combinations -----------

trimesters <- c("T01", "T02", "T03")

results_list <- list()

for (tr in trimesters) {
  sub <- mother %>% filter(corrected_visit_sample == tr)
  
  # Case/control -- T01 skipped explicitly (2 cases only, see header note)
  if (tr == "T01") {
    cat("\nSkipping case/control PLS-DA for T01 -- only 2 case samples.\n")
  } else {
    res <- run_plsda_permanova(sub, "casecat", "casecontrol", tr)
    if (!is.null(res)) results_list[[length(results_list) + 1]] <- res
  }
  
  # Maternal T1D -- run for all three trimesters
  res <- run_plsda_permanova(sub, "T1Dcat", "maternalT1D", tr)
  if (!is.null(res)) results_list[[length(results_list) + 1]] <- res
}

# ---- 4. Combined summary table -----------------------------------------------

summary_table <- bind_rows(results_list)
print(summary_table)

write_csv(summary_table, "multivariate_output/permanova_summary_all.csv")

cat("\nDone. Outputs saved in multivariate_output/:\n",
    "- plsda_indiv_<trimester>_<group>.jpeg     (PLS-DA ordination plot, supervised)\n",
    "- christmastree_<trimester>_<group>.jpeg   (top-20 loading/Christmas tree plot)\n",
    "- top20_<trimester>_<group>.csv            (top 20 contributing proteins, as a table)\n",
    "- pcoa_<trimester>_<group>.jpeg            (PCoA ordination plot, unsupervised)\n",
    "- permdisp_boxplot_<trimester>_<group>.jpeg (within-group dispersion boxplot)\n",
    "- permanova_summary_all.csv                (R2, PERMANOVA p, PCoA %var, PERMDISP p, across all runs)\n")