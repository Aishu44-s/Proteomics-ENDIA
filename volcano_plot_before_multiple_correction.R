# Run this once if ggrepel is not already installed:
# install.packages("ggrepel")

library(dplyr)
library(ggplot2)
library(ggrepel)

# ---------------------------------------------------------------------------
# Settings: update the input file path if needed
# ---------------------------------------------------------------------------

results_file <- paste0(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/",
  "NEW ANALYSIS/Helena_code/volcano_plot/",
  "volcano_case_control_by_trimester_no_maternal_T1D_data.csv"
)

out_dir <- paste0(
  "C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/",
  "NEW ANALYSIS/Helena_code/volcano_plot"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Change these three names only if your CSV uses different column names.
protein_col <- "Protein"
effect_col  <- "Effect"
raw_p_col   <- "PValue"

# ---------------------------------------------------------------------------
# Read and prepare results
# ---------------------------------------------------------------------------

results <- read.csv(results_file, check.names = FALSE)

required_columns <- c(protein_col, effect_col, raw_p_col)

if (!all(required_columns %in% names(results))) {
  stop(
    "These required columns were not found: ",
    paste(setdiff(required_columns, names(results)), collapse = ", "),
    "\n\nAvailable columns are:\n",
    paste(names(results), collapse = ", ")
  )
}

volcano_data <- results %>%
  transmute(
    Protein = as.character(.data[[protein_col]]),
    Effect = as.numeric(.data[[effect_col]]),
    RawPValue = as.numeric(.data[[raw_p_col]])
  ) %>%
  filter(
    is.finite(Effect),
    is.finite(RawPValue),
    RawPValue >= 0,
    RawPValue <= 1
  ) %>%
  mutate(
    RawPValue = pmax(RawPValue, .Machine$double.xmin),
    MinusLog10P = -log10(RawPValue),
    NominalAssociation = RawPValue < 0.10
  )

# Label only the ten proteins with the smallest raw p-values
label_data <- volcano_data %>%
  filter(NominalAssociation) %>%
  slice_min(RawPValue, n = 10, with_ties = FALSE)

# ---------------------------------------------------------------------------
# Exploratory volcano plot using RAW, unadjusted p-values
# ---------------------------------------------------------------------------

p_case_control_nominal <- ggplot(
  volcano_data,
  aes(x = Effect, y = MinusLog10P)
) +
  geom_point(
    aes(colour = NominalAssociation),
    size = 2.2,
    alpha = 0.80
  ) +
  geom_vline(
    xintercept = 0,
    colour = "grey50",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = -log10(0.10),
    colour = "red",
    linetype = "dotted",
    linewidth = 0.9
  ) +
  geom_text_repel(
    data = label_data,
    aes(label = Protein),
    size = 3.2,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    show.legend = FALSE
  ) +
  scale_colour_manual(
    values = c(
      "FALSE" = "grey70",
      "TRUE" = "#2C7FB8"
    ),
    labels = c(
      "FALSE" = "Raw p >= 0.10",
      "TRUE" = "Nominal raw p < 0.10"
    ),
    name = NULL
  ) +
  labs(
    title = "Exploratory case/control protein associations",
    subtitle = paste0(
      "No proteins met the FDR-adjusted q < 0.10 threshold"
    ),
    x = "Case/control effect estimate",
    y = expression(-log[10]("Raw p-value")),
    caption = paste0(
      "Red dotted line: nominal raw p = 0.10. ",
      "Results above this line were not taken forward because none ",
      "remained significant after multiple-testing correction."
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 11, lineheight = 1.1),
    plot.caption = element_text(hjust = 0, size = 9),
    legend.position = "top"
  )

print(p_case_control_nominal)

ggsave(
  filename = file.path(
    out_dir,
    "volcano_case_control_nominal_raw_p_lt_0.10.png"
  ),
  plot = p_case_control_nominal,
  width = 10,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)