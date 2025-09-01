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
d_initial <- readRDS(data_path) #%>%
  #amp_filter_samples(SampleSite == "Aalborg W")

# 2. Get unique sample sites from the metadata
# This identifies the different locations to be processed individually.
sample_sites <- unique(d_initial$metadata$SampleSite)

# 3. Create an empty list to store the final ampvis2 objects
# This list will hold the processed data for each site before exporting.
processed_ampvis_objects <- list()

# 4. Loop through each sample site to process them one by one
for (site in sample_sites) {
  # --- Filter by SampleSite ---
  cat(paste("Processing SampleSite:", site, "\n"))
  
  # Subset the ampvis2 object to only include data from the current site
  d_filtered <- amp_filter_samples(d_initial, SampleSite == site)
  # 
  # --- RCLR Transform ---
  # Apply the robust centered log-ratio transformation to the abundance data
  d_rclr <- rclr_transform(d_filtered)
  #
  # # --- Filter OTUs ---
  # As the data has been rclr transformed, we first filter the non-transformed data (`d_filtered`) as the rclr_subset cannot use filter_otus()
  # to get a list of OTUs to keep. Then we apply this filter to the transformed data (`d_rclr`) using the taxa_filter_vector.
  d_otu_filter <- filter_otus(d_filtered, filter_otus = 0.1)
  taxa_filter_vector <- as.vector(d_otu_filter$tax$OTU)

  # Subset the RCLR transformed data based on the taxa filter vector
  d_filtered_otus <- amp_subset_taxa(d_rclr, tax_vector = taxa_filter_vector)

  
  # Check if any data remains after filtering to avoid errors
  if (nrow(d_filtered_otus$abund) > 0 && ncol(d_filtered_otus$abund) > 0) {
    # Add the processed ampvis2 object to our list, naming it by the sample site
    processed_ampvis_objects[[site]] <- d_filtered_otus
  } else {
    # Print a message if a site has no data left after filtering
    cat(paste("  No data remaining for SampleSite:", site, "after filtering. Skipping.\n"))
  }
}

# 5. Export the processed abundance tables to a folder
output_folder <- here("output", "elsa_tables")
dir.create(output_folder, showWarnings = FALSE) # Create the folder if it doesn't exist

if (length(processed_ampvis_objects) > 0) {
  cat("\nExporting processed abundance tables...\n")
  
  # Loop through the list of processed ampvis2 objects
  for (site_name in names(processed_ampvis_objects)) {
    
    # Get the ampvis2 object for the current site
    amp_obj <- processed_ampvis_objects[[site_name]]
    
    # Create a clean filename by replacing spaces with underscores
    clean_site_name <- gsub(" ", "_", site_name)
    
    # Define the full path for the output file
    file_path <- here(output_folder, paste0(clean_site_name, "_rclr_abund"))
    
    # Export the abundance table
    amp_export_otutable(amp_obj, filename = file_path, sep = "\t", extension = "tsv", id = "SampleDate")
    
    cat(paste("  Saved table for", site_name, "to:", file_path, "\n"))
  }
  cat("\nExport complete.\n")
} else {
  cat("\nNo data was processed. No tables were exported.\n")
}

# 6. Load, reformat, sort, and re-save the tables
cat("\nRe-formatting and sorting tables for elsa...\n")
sorted_output_folder <- here("output", "elsa_tables") 
#dir.create(sorted_output_folder, showWarnings = FALSE, recursive = TRUE)

# Get a list of the files that were just created
exported_files <- list.files(output_folder, pattern = "\\.tsv$", full.names = TRUE)

for (file_path in exported_files) {
  # Read the exported table
  # `check.names = FALSE` prevents R from altering the date format in the column names
  abund_table <- read.delim(file_path, sep = "\t", row.names = 1, check.names = FALSE)%>%
    select(!(Kingdom:Species))
  
  # Convert column names to Date objects
  date_colnames <- colnames(abund_table)
  dates <- as.Date(date_colnames)
  
  # Sort the table columns in ascending order based on the dates
  sorted_table <- abund_table[, order(dates)]
  
  # Rename columns to the format tXr1 (e.g., t1r1, t2r1, ...)
  # This assumes each column is a unique timepoint with one replicate.
  colnames(sorted_table) <- paste0("t", 1:ncol(sorted_table), "r1")
  
  # Define the new file path for the sorted table
  new_file_name <- basename(file_path)
  new_file_path <- here(sorted_output_folder, paste0("sorted_", new_file_name))
  
  # To ensure the output format is correct for elsa, we'll add the OTU IDs as the first column
  # and then write the table without row names.
  # The first column header will be "#OTU_ID" and the rest will be the new tXr1 names.
  output_df <- data.frame("#OTU_ID" = rownames(sorted_table), sorted_table, check.names = FALSE)
  
  # Write the sorted table to the new folder, using 'na' for missing values
  write.table(
    output_df, 
    file = new_file_path, 
    sep = "\t", 
    quote = FALSE, 
    row.names = FALSE,
    na = "na"
  )
  
  cat(paste("  Saved formatted table to:", new_file_path, "\n"))
}

cat("\nTable re-formatting complete.\n")
