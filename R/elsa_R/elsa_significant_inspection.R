# Load necessary libraries
# Make sure you have these installed: install.packages(c("here", "readr", "dplyr", "purrr"))
library(here)
library(readr)
library(dplyr)
library(purrr)

# 1. Define the path to your eLSA result files
results_path <- here("output", "elsa_results")

# 2. Get a list of all files ending with .lsa
all_lsa_files <- list.files(
  path = results_path,
  pattern = "\\.lsa$", # This gets ALL files ending in .lsa
  full.names = TRUE
)

norm_off <- FALSE # TRUE = Load norm_off.lsa files

# --- NEW: Exclude the norm_off files ---
# From the full list, remove any files that end with "norm_off.lsa"
file_list <- grep("norm_off\\.lsa$", all_lsa_files, value = TRUE, invert = ifelse(norm_off == FALSE, TRUE, FALSE))

# 3. Loop over the filtered file list and read each one
# This part remains the same; it will now operate on the corrected file list.
all_lsa_data <- map_dfr(file_list, ~{
  # Extract the sample site name from the file path
  sample_site <- basename(.x) %>%
    gsub(ifelse(norm_off == FALSE, "\\_rclr_abund.lsa$", "\\_rclr_abund_norm_off.lsa$"), "", .) # A more general way to clean the name
  
  # Read the .lsa file
  read_tsv(.x, col_types = cols(.default = "c")) %>%
    mutate(
      SampleSite = sample_site,
      .before = 1
    ) %>%
    filter(Delay != 0, 
           Q <= 0.01)
})

# 4. Safely convert columns to numeric (error-proof method)
expected_numeric_cols <- c("LS", "Delay", "Q", "P", "PCC", "SCC", "SPCC", "SSCC")
cols_to_convert <- intersect(expected_numeric_cols, colnames(all_lsa_data))

if (length(cols_to_convert) > 0) {
  all_lsa_data <- all_lsa_data %>%
    mutate(across(
      all_of(cols_to_convert),
      as.numeric
    ))
}

# 5. Final inspection
print(paste("Successfully combined", length(file_list), "files (after excluding 'norm_off' files)."))
print(paste("The final data frame has", nrow(all_lsa_data), "rows and", ncol(all_lsa_data), "columns."))
print("Columns converted to numeric:")
print(cols_to_convert)

# Display the first few rows of the result
sig_interactions_df <- all_lsa_data %>%
  group_by(SampleSite) %>%
  summarise(interactions = n()) %>%
  mutate(sum = sum(interactions))

head(sig_interactions_df)

# Optional: Save the combined data frame to a file
# write_csv(all_lsa_data, here("output", "combined_lsa_results_filtered.csv"))