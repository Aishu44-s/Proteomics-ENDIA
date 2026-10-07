library(ggplot2)
library(dplyr)

# File locations
project_dir <- "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics"
analysis_dir <- file.path(project_dir, "NEW ANALYSIS")
merge_dir <- file.path(analysis_dir, "Merging data NEW")
filter_dir <- file.path(analysis_dir, "Missing values and filtering NEW")

# Read the datasets and protein lists
merged <- read.csv(
  file.path(merge_dir, "merged_meta_protein.csv"),
  check.names = FALSE
)

filtered <- read.csv(
  file.path(filter_dir, "protein_filtered_80pct.csv"),
  check.names = FALSE
)

final_data <- read.csv(
  file.path(filter_dir, "protein_filtered_80pct_outlier_removed.csv"),
  check.names = FALSE
)

raw_proteins <- readRDS(file.path(project_dir, "protein_column_names.rds"))
filtered_proteins <- readRDS(
  file.path(filter_dir, "protein_column_names_post_filter.rds")
)
outlier_ids <- readRDS(file.path(filter_dir, "outlier_removed_sample_ids.rds"))

# Calculate sample and protein counts
n_metadata_rows <- nrow(merged)
n_unique_samples <- n_distinct(merged$SampleID.2)
n_samples_after_filter <- nrow(filtered)
n_outliers <- length(outlier_ids)
n_final_samples <- nrow(final_data)

n_raw_proteins <- length(raw_proteins)
n_retained_proteins <- length(filtered_proteins)
n_not_retained <- length(setdiff(raw_proteins, filtered_proteins))

cat("Metadata-linked rows:", n_metadata_rows, "\n")
cat("Unique physical samples:", n_unique_samples, "\n")
cat("Outlier samples removed:", n_outliers, "\n")
cat("Final samples:", n_final_samples, "\n")
cat("Proteins in raw dataset:", n_raw_proteins, "\n")
cat("Proteins not in filtered list:", n_not_retained, "\n")
cat("Proteins in filtered list:", n_retained_proteins, "\n")

# Draw a simple vertical flow diagram
make_flow_plot <- function(labels, title, output_file) {
  steps <- data.frame(
    label = labels,
    x = 1,
    y = rev(seq_along(labels))
  )

  arrows <- data.frame(
    x = 1,
    xend = 1,
    y = steps$y[-nrow(steps)] - 0.35,
    yend = steps$y[-1] + 0.35
  )

  p <- ggplot() +
    geom_segment(
      data = arrows,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = grid::arrow(length = grid::unit(0.18, "inches")),
      linewidth = 0.7,
      colour = "#666666"
    ) +
    geom_label(
      data = steps,
      aes(x = x, y = y, label = label),
      fill = "#F5F1E8",
      colour = "#333333",
      label.size = 0.5,
      label.padding = grid::unit(0.28, "lines"),
      size = 4.5,
      lineheight = 1.05
    ) +
    coord_cartesian(xlim = c(0.4, 1.6), ylim = c(0.5, nrow(steps) + 0.5)) +
    labs(title = title) +
    theme_void(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        colour = "#222222",
        margin = margin(b = 15)
      )
    )

  ggsave(output_file, p, width = 8, height = 2.1 * nrow(steps), dpi = 600)
  p
}

# Sample flow
sample_flow <- make_flow_plot(
  labels = c(
    paste0("Metadata-linked rows\nn = ", n_metadata_rows,
           "\nIncludes repeated nested-control records"),
    paste0("Unique physical maternal stool samples\nn = ", n_unique_samples),
    paste0("Outlier assessment\n", n_outliers, " samples excluded"),
    paste0("Samples entering imputation and downstream analysis\nn = ",
           n_final_samples)
  ),
  title = "Sample flow",
  output_file = "sample_flow_diagram.png"
)

# Protein flow
raw_proteins <- readRDS(
  file.path(project_dir, "protein_column_names.rds")
)

filtered_proteins <- readRDS(
  file.path(filter_dir, "protein_column_names_post_filter.rds")
)

n_raw <- length(raw_proteins)
n_retained <- length(filtered_proteins)
n_not_retained <- n_raw - n_retained

cat("Raw protein count:", n_raw, "\n")
cat("Post-filter protein count:", n_retained, "\n")
cat("Count difference:", n_not_retained, "\n")
cat("Protein names shared between lists:",
    length(intersect(raw_proteins, filtered_proteins)), "\n")

stopifnot(n_not_retained >= 0)

nodes <- data.frame(
  x = c(2, 1, 3),
  y = c(3.4, 1.5, 1.5),
  label = c(
    paste0("Proteins quantified in raw dataset\nn = ", n_raw),
    paste0("Not retained after filtering/QC\nn = ", n_not_retained),
    paste0("Retained in post-filter dataset\nn = ", n_retained)
  )
)

arrows <- data.frame(
  x = c(1.9, 2.1),
  y = c(3.0, 3.0),
  xend = c(1.15, 2.85),
  yend = c(1.95, 1.95)
)

protein_flow <- ggplot() +
  geom_segment(
    data = arrows,
    aes(x = x, y = y, xend = xend, yend = yend),
    arrow = grid::arrow(length = grid::unit(0.18, "inches")),
    linewidth = 0.8,
    colour = "#555555"
  ) +
  geom_label(
    data = nodes,
    aes(x = x, y = y, label = label),
    fill = "#F5F1E8",
    colour = "#333333",
    label.size = 0.5,
    size = 5,
    lineheight = 1.1,
    label.padding = grid::unit(0.35, "lines")
  ) +
  coord_cartesian(xlim = c(0.3, 3.7), ylim = c(0.8, 4.1)) +
  labs(title = "Protein filtering flow") +
  theme_void(base_size = 14) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      colour = "#222222",
      margin = margin(b = 18)
    ),
    plot.background = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

ggsave(
  file.path(analysis_dir, "workflow", "protein_flow_diagram.png"),
  protein_flow,
  width = 10,
  height = 6,
  dpi = 600,
  bg = "white"
)