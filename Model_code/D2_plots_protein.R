library(ggplot2)

# ── Update resp to all 451 proteins ──────────────────────────────
resp  <- names(mother)[grepl("^(sp|tr)\\.", names(mother))]
resp3 <- resp   # use column names as display labels

# ── Create output folder ──────────────────────────────────────────
dir.create("violin_plots", showWarnings = FALSE)

cat("Starting full run for", length(resp), "proteins...\n")

for(i in 1:length(resp)){
  
  tmp <- mother
  w   <- which(names(tmp) == resp[i])
  names(tmp)[w] <- "btmean"
  tmp <- tmp[!is.na(tmp$btmean), ]
  
  tmp$casecat <- factor(tmp$casecat, levels = c(0,1),
                        labels = c("Control","Case"))
  
  tmp$caset1dcat <- as.character(tmp$caset1dcat)
  tmp$caset1dcat <- factor(tmp$caset1dcat,
                           levels = c("IA_0_T1D_0","IA_0_T1D_1","IA_1_T1D_0","IA_1_T1D_1"),
                           labels = c("Control Without T1D","Control With T1D",
                                      "Case Without T1D","Case With T1D"))
  
  g <- ggplot(data = tmp, aes(x = caset1dcat, y = btmean, fill = caset1dcat)) +
    geom_violin(aes(group = caset1dcat)) +
    geom_point(position = position_jitter(seed = 1, width = 0.2)) +
    stat_summary(fun = median, geom = "point", size = 4, color = "purple") +
    ylab(resp3[i]) + xlab("IA status") +
    theme_classic() + theme(legend.position = "none")
  
  g <- g + scale_fill_manual(values = c("light blue","lightblue3","pink","pink3"))
  
  g1 <- g + labs(color = "Controls") +
    theme_classic() + guides(fill = "none")
  
  # Save with safe filename
  safe_name <- gsub("[^A-Za-z0-9_\\-]", "_", resp3[i])
  ggsave(paste0("violin_plots/", safe_name, "_t1d_graphviolin.jpeg"),
         plot = g1, width = 8, height = 6, dpi = 150)
  
  # Progress tracker — prints every 50 proteins
  if(i %% 50 == 0) cat("Progress:", i, "/", length(resp), "done\n")
}

cat("ALL DONE —", length(resp), "plots saved to violin_plots/ folder\n")
