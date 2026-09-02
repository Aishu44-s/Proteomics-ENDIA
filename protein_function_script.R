library(dplyr)
library(readr)
library(httr)
library(stringr)

# ---- 1. Load your hit list --------------------------------------------------
# Point this at whichever protein list you want annotated -- e.g. the overlap
# file from the q-value script, or the full qvalue_long file filtered to hits.

hits <- read_csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Helena_code/qvalue_plots/visit_t1dvisit_overlap_proteins.csv", show_col_types = FALSE)
# hits has columns: protein, category

# ---- 2. Extract UniProt accessions from the protein column names -----------
# Protein names look like "sp.P01023.A2MG" -> accession is the middle part.

hits <- hits %>%
  mutate(accession = str_split_fixed(protein, "\\.", 3)[, 2])

accessions <- unique(hits$accession)
cat("Number of unique accessions to look up:", length(accessions), "\n")

# ---- 3. Query UniProt REST API in batches -----------------------------------
# UniProt's search endpoint accepts a boolean OR query across accessions and
# returns a TSV. Batching avoids overly long URLs (UniProt caps ~ a few
# hundred accessions comfortably per request; 100 is a safe chunk size).

batch_size <- 100
batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))

fetch_batch <- function(acc_batch) {
  query <- paste(sprintf("accession:%s", acc_batch), collapse = " OR ")
  url <- "https://rest.uniprot.org/uniprotkb/search"
  
  resp <- GET(url, query = list(
    query  = query,
    fields = "accession,gene_names,protein_name,cc_function,organism_name",
    format = "tsv",
    size   = 500
  ))
  
  if (status_code(resp) != 200) {
    warning("UniProt request failed with status ", status_code(resp))
    return(NULL)
  }
  
  read_tsv(content(resp, as = "text", encoding = "UTF-8"), show_col_types = FALSE)
}

annotation_list <- lapply(batches, fetch_batch)
annotations <- bind_rows(annotation_list)

names(annotations) <- c("accession", "gene_names", "protein_name", "function_description", "organism")

# ---- 4. Merge annotations back onto your hit list ---------------------------

hits_annotated <- hits %>%
  left_join(annotations, by = "accession")

write_csv(hits_annotated, "protein_function_annotations.csv")

cat("\nDone. Annotated table written to protein_function_annotations.csv\n")
cat("Columns: protein, category, accession, gene_names, protein_name, function_description, organism\n")

# ---- 5. (Optional) merge with q-values too, if you want function + q-value
#         side by side in one table -----------------------------------------
# qlong <- read_csv("qvalue_long_form_adj_int_rm.csv", show_col_types = FALSE)
# qlong <- qlong %>% mutate(accession = str_split_fixed(protein, "\\.", 3)[, 2])
# full_table <- qlong %>% left_join(annotations, by = "accession")
# write_csv(full_table, "qvalue_with_function_annotations.csv")