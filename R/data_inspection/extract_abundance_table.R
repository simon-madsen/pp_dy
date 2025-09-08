# Load necessary libraries
library(ampvis2)
library(here) # For robust path management
library(readr) # For efficient file reading and writing

# Load the dataset
d <- readRDS(here("data", "d_initial_Simon_subset.RDS"))

# --- 1. Set up the output directory ---
output_dir <- here("data", "abundance_tables")

# Create the directory if it doesn't already exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Created directory:", output_dir, "\n")
}

# --- 2. Loop, Export, Process, and Resave ---
sample_sites <- unique(d$metadata$SampleSite)
processed_tables <- list()

for (site in sample_sites) {
  cat("Processing site:", site, "...\n")
  
  # Subset the ampvis2 object for the current site
  d_site <- amp_filter_samples(d, SampleSite == site, normalise = FALSE)
  
  # Define the full path for the output file
  file_path <- file.path(output_dir, paste0(site, "_abundance_table.tsv"))
  
  # a) Export the initial abundance table
  amp_export_otutable(
    d_site,
    id = "SampleDate",
    filename = file_path,
    extension = ""
  )
  
  # b) Re-open the file
  ab_table <- read_delim(file_path, delim = "\t", show_col_types = FALSE)
  
  # Check for sufficient data and skip if necessary
  if (ncol(ab_table) <= 1) {
    cat("  -> Skipped. No sample data to process.\n")
    next
  }
  
  # c) Modify the reopened table
  if (ncol(ab_table) == 2) {
    # Handle the simple case of a single sample date
    modified_table <- ab_table
    colnames(modified_table)[2] <- 0
  } else {
    # Handle the case with multiple sample dates
    date_colnames <- colnames(ab_table)[-1]
    dates <- as.Date(date_colnames)
    sorted_indices <- order(dates)
    
    modified_table <- ab_table[, c(1, sorted_indices + 1)]
    
    sorted_dates_str <- colnames(modified_table)[-1]
    sorted_dates <- as.Date(sorted_dates_str)
    numeric_days <- as.numeric(sorted_dates) - as.numeric(sorted_dates[1])
    
    colnames(modified_table) <- c("OTU", numeric_days)
  }
  
  # d) Resave the MODIFIED table to the SAME path, overwriting the original
  write_tsv(modified_table, file_path)
  
  # e) Store the final table in the R session list for immediate use
  processed_tables[[site]] <- modified_table
  
  cat("  -> Done. Modified table saved to:", file_path, "\n")
}

# --- 3. Display an Example ---
# The list `processed_tables` contains the final data for your current R session.
# The files in `data/abundance_tables` now contain the versions with numeric day columns.
if (length(processed_tables) > 0) {
  cat("\n--- Example Output for the first processed site ---\n")
  dplyr::glimpse(processed_tables[[1]])
}