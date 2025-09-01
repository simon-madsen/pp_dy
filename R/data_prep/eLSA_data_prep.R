# Load necessary libraries
# If you don't have these, install them with: install.packages(c("tidyverse", "ampvis2"))
library(tidyverse)
library(ampvis2)
library(here)
# Source the custom RCLR transformation function
# Make sure the 'rclr_transform.R' file is in your working directory.
# If it's in a different directory, provide the full path to the file.
# For example: source("C:/Users/Simon/Documents/Msc project/msc_R/rclr_transform.R")
source(here("R", "functions", "rclr_transform.R"))

# 1. Load your ampvis2 object from the specified .rds file
# Ensure the file path is correct and accessible from your R environment.
# Using double backslashes (\\) or a single forward slash (/) is recommended in R for file paths.
data_path <- here("data", "d_initial_Simon_subset.rds")
d_initial <- readRDS(data_path) %>%
  amp_filter_samples(SampleSite == "Aalborg W")

# 2. Get unique sample sites from the metadata
# This identifies the different locations to be processed individually.
sample_sites <- unique(d_initial$metadata$SampleSite)

# 3. Define the output folder and create it if it doesn't exist
# The `recursive = TRUE` argument ensures that both 'output' and 'elsa_tables' are created if needed.
output_folder <- here("output", "elsa_tables")
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)


# 4. Loop through each sample site to process and export them one by one
for (site in sample_sites) {
  cat(paste("Processing SampleSite:", site, "\n"))
  
  # --- Step A: Filter and Transform Data ---
  d_filtered <- amp_filter_samples(d_initial, SampleSite == site)
  d_rclr <- rclr_transform(d_filtered)
  
  # Filter non-transformed data to get a list of ASVs to keep
  d_otu_filter <- filter_otus(d_filtered, filter_otus = 0.1)
  taxa_filter_vector <- as.vector(d_otu_filter$tax$OTU)
  
  # Subset the RCLR transformed data
  d_processed <- amp_subset_taxa(d_rclr, tax_vector = taxa_filter_vector)
  
  # --- Step B: Check for Data and Extract Abundance Table ---
  if (nrow(d_processed$abund) > 0 && ncol(d_processed$abund) > 0) {
    abund_table <- d_processed$abund
    
    # --- Step C: Sort by Date ---
    # Create a lookup for sample IDs and their dates from the metadata
    meta_subset <- d_processed$metadata
    date_lookup <- setNames(as.Date(meta_subset$SampleDate), meta_subset$SampleID)
    
    # Get the dates for the columns in the abundance table
    col_dates <- date_lookup[colnames(abund_table)]
    
    # Sort the table columns by date in ascending order
    sorted_table <- abund_table[, order(col_dates)]
    
    # --- Step D: Format for eLSA ---
    # Rename columns to the format tXr1
    colnames(sorted_table) <- paste0("t", 1:ncol(sorted_table), "r1")
    
    # Create the final data frame with the required "#OTU_ID" header
    output_df <- data.frame("#OTU_ID" = rownames(sorted_table), sorted_table, check.names = FALSE)
    
    # --- Step E: Export the Final Table ---
    # Define the final output file path
    clean_site_name <- gsub(" ", "_", site)
    file_path <- here(output_folder, paste0(clean_site_name, "_rclr_abund.tsv"))
    
    # Write the final table
    write.table(
      output_df, 
      file = file_path, 
      sep = "\t", 
      quote = FALSE, 
      row.names = FALSE,
      na = "na"
    )
    cat(paste("  Successfully saved formatted table to:", file_path, "\n"))
    
  } else {
    cat(paste("  No data remaining for SampleSite:", site, "after filtering. Skipping.\n"))
  }
}

cat("\nScript finished. All files are in the '", output_folder, "' directory.\n")


