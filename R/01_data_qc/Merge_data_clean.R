# Merge ENDIA metadata with maternal stool protein data

library(dplyr)

# Load input files
meta <- read.csv(
  "C:/Users/Dell/Documents/BIOX7011/Data/Filtered_maternal_metadata.csv",
  check.names = FALSE
)

raw <- read.csv(
  "C:/Users/Dell/Documents/BIOX7011/Data/raw_maternal_human_protein.csv",
  check.names = FALSE
)

names(raw)[1] <- "SampleID_raw"

# Save the protein column names for downstream scripts
protein_col_names <- setdiff(names(raw), c("SampleID_raw", "Plate"))
saveRDS(protein_col_names, "protein_column_names.rds")

cat("Metadata rows:", nrow(meta), "\n")
cat("Raw protein rows:", nrow(raw), "\n")
cat("Protein columns:", length(protein_col_names), "\n")
cat("Reused-control metadata rows:", sum(grepl("_a$", meta$SampleID)), "\n")
cat("Raw IDs ending in _a:", sum(grepl("_a$", raw$SampleID_raw)), "\n")

# Check that raw sample IDs are unique and metadata samples have raw data
stopifnot(!anyDuplicated(raw$SampleID_raw))

missing_ids <- setdiff(meta$SampleID.2, raw$SampleID_raw)
if (length(missing_ids) > 0) {
  stop("These metadata SampleID.2 values have no raw-data match: ",
       paste(missing_ids, collapse = ", "))
}

# Join on the physical sample ID. Reused controls are copied to each matching row.
merged <- left_join(meta, raw, by = c("SampleID.2" = "SampleID_raw"))

stopifnot(nrow(merged) == nrow(meta))
stopifnot(all(protein_col_names %in% names(merged)))

if ("Plate.y" %in% names(merged) && anyNA(merged$Plate.y)) {
  stop("Some metadata rows did not match a raw protein sample.")
}

# Save the full dataset for modelling and one-row-per-sample data for QC
write.csv(merged, "merged_meta_protein.csv", row.names = FALSE)

unique_samples <- distinct(merged, SampleID.2, .keep_all = TRUE)
write.csv(unique_samples, "unique_samples_meta_protein.csv", row.names = FALSE)

cat("Merged rows:", nrow(merged), "\n")
cat("Unique physical samples:", nrow(unique_samples), "\n")
cat("Saved merged_meta_protein.csv, unique_samples_meta_protein.csv, and protein_column_names.rds\n")
