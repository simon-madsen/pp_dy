# Load necessary libraries
# If you don't have these, install them with: install.packages(c("tidyverse", "ampvis2"))
library(tidyverse)
library(ampvis2)

#' Robust Centered Log-Ratio (RCLR) Transformation
#'
#' This function applies a robust centered log-ratio transformation to the 
#' abundance table of an ampvis2 object. The transformation is applied per-sample (column-wise).
#'
#' @param ampvis_object An object of class 'ampvis2'.
#'
#' @return A new 'ampvis2' object with the abundance data transformed, preserving the original structure.
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data(MiDAS)
#'
#' # Apply the rclr transformation
#' MiDAS_rclr <- rclr_transform(MiDAS)
#'
#' # Check the transformed data
#' head(MiDAS_rclr$abund[, 1:5])
#' }

rclr_transform <- function(ampvis_object) {
  # Ensure the input is a valid ampvis2 object
  if (!inherits(ampvis_object, "ampvis2")) {
    stop("Error: Input must be a valid ampvis2 object.")
  }
  
  # >> FIX: Preserve original data frame structure (rownames and colnames) <<
  # Store original row and column names to re-apply them later.
  original_rownames <- rownames(ampvis_object$abund)
  original_colnames <- colnames(ampvis_object$abund)
  
  # Extract only the numeric abundance data for transformation.
  # This assumes OTU identifiers are in rownames.
  abund_table <- ampvis_object$abund %>%
    select_if(is.numeric)
  
  # Define the rclr transformation function to be applied to each column (sample)
  rclr_sample <- function(sample_counts) {
    # Keep only the observed (non-zero) counts for the geometric mean calculation
    observed_counts <- sample_counts[sample_counts > 0]
    
    # If a sample has no observed counts, return it as is (all zeros)
    if (length(observed_counts) == 0) {
      return(sample_counts)
    }
    
    # Calculate the geometric mean of the observed counts
    gmean <- exp(mean(log(observed_counts)))
    
    # Apply the log-ratio transformation: log(count / geometric_mean)
    # For zero counts, the result will be -Inf. We will replace these with NA.
    transformed_counts <- log(sample_counts / gmean)
    transformed_counts[is.infinite(transformed_counts)] <- NA
    
    return(transformed_counts)
  }
  
  # Apply the rclr_sample function to each column (sample) of the abundance table
  transformed_abund <- as.data.frame(lapply(abund_table, rclr_sample))
  
  # Replace NA values with 0 to prevent errors in downstream ampvis2 functions
  transformed_abund[is.na(transformed_abund)] <- 0
  
  # >> FIX: Restore original row and column names to prevent modification <<
  rownames(transformed_abund) <- original_rownames
  colnames(transformed_abund) <- original_colnames
  
  # Create a new ampvis2 object to store the results
  ampvis_transformed <- ampvis_object
  
  # Replace the original abundance data with the transformed data
  ampvis_transformed$abund <- transformed_abund
  
  # Return the new ampvis2 object with the transformed data
  return(ampvis_transformed)
}

# --- Example Usage ---
# 
# # # 1. Load some example data from the ampvis2 package
# data(MiDAS)
# #
# taxa_filter <- MiDAS %>%
#   filter_otus(filter_otus = 1)
# 
# taxa_filter_vector <- as.vector(taxa_filter$tax$OTU)
# 
# 
# # # 2. Apply the rclr transformation to the MiDAS dataset
# MiDAS_rclr <- rclr_transform(MiDAS)
# 
# MiDAS_rclr_subset <- amp_subset_taxa(MiDAS_rclr, tax_vector = taxa_filter_vector)
# 
# #
# # # 3. View the first few rows and columns of the original abundance table
# cat("Original Abundance Table:\n")
# print(head(MiDAS$abund[, 1:6]))
# #
# # # 4. View the first few rows and columns of the transformed abundance table
# # # Note the NA values where original counts were 0.
# cat("\nTransformed (RCLR) Abundance Table:\n")
# print(head(MiDAS_rclr_subset$abund[, 1:6]))

