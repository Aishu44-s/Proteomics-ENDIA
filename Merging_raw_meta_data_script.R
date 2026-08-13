#===================================================================

library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(limma)

#-------------------------------------------------------------------
# 1. Load data 
#-------------------------------------------------------------------

meta <- read.csv("C:/Users/Dell/Documents/BIOX7011/Data/Filtered_maternal_metadata.csv", header = TRUE)
raw  <- read.csv("C:/Users/Dell/Documents/BIOX7011/Data/raw_maternal_human_protein.csv",  header = TRUE)

# The raw file's first column (sample name) imports with a blank/auto
# header — rename it explicitly so it's unambiguous downstream
names(raw)[1] <- "SampleID_raw"

cat("Metadata rows:", nrow(meta), "\n")
cat("Raw protein rows:", nrow(raw), "\n")

#-------------------------------------------------------------------
# 1b. *** NEW *** Capture the definitive protein column list here,
#     while we still have unambiguous ground truth: everything in
#     `raw` except the sample ID and the Plate column is a protein.
#-------------------------------------------------------------------

protein_col_names <- setdiff(names(raw), c("SampleID_raw", "Plate"))

cat("\nProtein columns identified in raw file:", length(protein_col_names), "\n")
cat("  sp| (Swiss-Prot):", sum(str_starts(protein_col_names, "sp\\.")), "\n")
cat("  tr| (TrEMBL):    ", sum(str_starts(protein_col_names, "tr\\.")), "\n")
cat("  other/unexpected:", sum(!str_starts(protein_col_names, "sp\\.") &
                                 !str_starts(protein_col_names, "tr\\.")), "\n")
# If "other/unexpected" is not 0, print them so you can eyeball what they are
# before trusting the pipeline downstream
unexpected <- protein_col_names[!str_starts(protein_col_names, "sp\\.") &
                                  !str_starts(protein_col_names, "tr\\.")]
if (length(unexpected) > 0) {
  cat("  Unexpected column names (inspect these):\n")
  print(unexpected)
}

saveRDS(protein_col_names, "protein_column_names.rds")
cat("Saved: protein_column_names.rds (source of truth for all downstream scripts)\n")

#-------------------------------------------------------------------
# 2. Sanity checks before merging
#-------------------------------------------------------------------

# How many metadata rows are "_a" (reused-control) rows?
n_a <- sum(str_detect(meta$SampleID, "_a$"))
cat("\nMetadata rows with '_a' suffix:", n_a, "\n")

# Raw file should have NO "_a" names at all (physical samples only)
n_a_raw <- sum(str_detect(raw$SampleID_raw, "_a$"))
cat("Raw rows with '_a' suffix (should be 0):", n_a_raw, "\n")

# Every metadata SampleID.2 should exist in the raw file
missing_in_raw <- setdiff(meta$SampleID.2, raw$SampleID_raw)
cat("meta$SampleID.2 values with NO match in raw data:", length(missing_in_raw), "\n")
if (length(missing_in_raw) > 0) print(missing_in_raw)

# Raw sample names should be unique (one row per physical sample)
cat("Duplicate names in raw data (should be 0):",
    sum(duplicated(raw$SampleID_raw)), "\n")

#-------------------------------------------------------------------
# 3. Merge — join on SampleID.2 (plain name), which exists in raw
#-------------------------------------------------------------------

merged <- meta %>%
  left_join(raw, by = c("SampleID.2" = "SampleID_raw"))

cat("\nRows in merged data:", nrow(merged), "\n")
cat("Rows in metadata:   ", nrow(meta), "\n")

# These should match exactly. If merged has MORE rows than meta,
# you have an accidental many-to-many match somewhere upstream.
stopifnot(nrow(merged) == nrow(meta))

# No metadata row should have failed to find a protein match
# (Plate is a raw-data column that will be NA only if the join failed)
n_unmatched <- sum(is.na(merged$Plate.y))
cat("Metadata rows with NO matching protein data:", n_unmatched, "\n")
if (n_unmatched > 0) {
  print(merged$SampleID[is.na(merged$Plate.y)])
}

# *** NEW *** Confirm every expected protein column actually made it
# into `merged` intact (protects against silent column-name mangling
# differences between how `raw` was read here vs. how it might be
# re-read elsewhere)
missing_prot_cols <- setdiff(protein_col_names, names(merged))
cat("\nExpected protein columns missing from merged data (should be 0):",
    length(missing_prot_cols), "\n")
if (length(missing_prot_cols) > 0) print(missing_prot_cols)

# Confirm the '_a' rows really did get duplicated protein values
# (spot check: one 'a' row and its base row should be identical)
check_id <- meta$SampleID[str_detect(meta$SampleID, "_a$")][1]
base_id  <- str_remove(check_id, "_a$")

a_vals    <- merged %>% filter(SampleID == check_id)    %>% select(all_of(protein_col_names))
base_vals <- merged %>% filter(SampleID == base_id)     %>% select(all_of(protein_col_names))

cat("\nSpot check —", check_id, "vs", base_id, "\n")
cat("Protein values identical:", isTRUE(all.equal(a_vals, base_vals, check.attributes = FALSE)), "\n")

#-------------------------------------------------------------------
# 4. Save the merged dataset
#-------------------------------------------------------------------

write.csv(merged, "merged_meta_protein.csv", row.names = FALSE)
cat("\nSaved: merged_meta_protein.csv (", nrow(merged), "rows x", ncol(merged), "cols )\n")

#-------------------------------------------------------------------
# 5. A note for the QC / missingness / imputation / normalisation
#    steps that come next
#-------------------------------------------------------------------
# unique_samples does not contain the "_a" samples as it will affect imputation and normalisation.
# However it will be used for the model

unique_samples <- merged %>%
  distinct(SampleID.2, .keep_all = TRUE)

cat("Unique physical samples for QC/imputation/normalisation:",
    nrow(unique_samples), "\n")
