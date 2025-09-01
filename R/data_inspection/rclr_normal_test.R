# Load necessary libraries
# install.packages("tidyverse")
library(tidyverse)
library(here)
# The folder containing the sorted and formatted tables for LSA.
# This should be the final output folder from the 'lsa_analysis_script'.
input_folder <- here("data", "elsa_tables")

# The significance level (alpha) for the Shapiro-Wilk test.
# A p-value above this level suggests the data is normally distributed.
alpha <- 0.05

# --- Script Execution ---

# Check if the input directory exists
if (!dir.exists(input_folder)) {
  stop(paste("Error: Input directory '", input_folder, "' not found. Please run the data preparation script first.", sep=""))
}

# Get a list of all the .tsv files in the input folder
# The pattern now looks for the 'sorted_' prefix added in the previous script.
files_to_test <- list.files(input_folder, pattern = "^sorted_.*\\.tsv$", full.names = TRUE)

if (length(files_to_test) == 0) {
  stop(paste("No sorted .tsv files found in '", input_folder, "'. Check the file naming and location.", sep=""))
}

cat(paste("Found", length(files_to_test), "files to test for normality.\n\n"))

# Loop through each file to perform the normality tests
for (file_path in files_to_test) {
  cat(paste("--- Testing File:", basename(file_path), "---\n"))
  
  # Read the data table. `check.names = FALSE` is important to keep column names as is.
  # The first column is the OTU identifier.
  data <- read.delim(file_path, sep = "\t", check.names = FALSE, colClasses = c("character", rep("numeric", ncol(read.delim(file_path, nrows = 1)) - 1)))
  
  # Initialize counters for the results
  normal_count <- 0
  not_normal_count <- 0
  
  # Check if there are numeric columns to test
  if (ncol(data) < 2) {
    cat("  No data columns found to test in this file. Skipping.\n\n")
    next # Skip to the next file
  }
  
  # Loop through each numeric column (from the second column onwards)
  for (col_name in colnames(data)[-1]) {
    # Extract the column as a numeric vector
    column_data <- data[[col_name]]
    
    # Remove any 'na' values before testing, as shapiro.test requires at least 3 non-missing values
    #column_data <- na.omit(column_data)
    
    # Shapiro-Wilk test requires at least 3 data points
    if (length(column_data) >= 3) {
      # Perform the Shapiro-Wilk test for normality
      test_result <- shapiro.test(column_data)
      
      # Check the p-value against our significance level
      if (test_result$p.value > alpha) {
        # If p > alpha, we fail to reject the null hypothesis (data is normally distributed)
        normal_count <- normal_count + 1
      } else {
        # If p <= alpha, we reject the null hypothesis (data is not normally distributed)
        not_normal_count <- not_normal_count + 1
      }
    } else {
      # If there are not enough data points, we can't test it.
      cat(paste("  Column", col_name, "has fewer than 3 non-NA values. Skipping test for this column.\n"))
    }
  }
  
  # Print the summary for the current file
  cat(paste("  Total columns tested:", normal_count + not_normal_count, "\n"))
  cat(paste("  Columns assumed to be normally distributed (p >", alpha, "):", normal_count, "\n"))
  cat(paste("  Columns NOT assumed to be normally distributed (p <=", alpha, "):", not_normal_count, "\n\n"))
}

cat("--- Normality testing complete. ---\n")
