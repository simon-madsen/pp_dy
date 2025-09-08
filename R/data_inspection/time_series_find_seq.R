# --- --- --- --- --- ---
#      1. Load Libraries
# --- --- --- --- --- ---
library(ampvis2)
library(here)


# --- --- --- --- --- ---
#      2. Function Definitions
# --- --- --- --- --- ---

#' @title Find Longest, Most Flexible Consecutive Sequence (BLAST-like)
#' @description This function identifies the longest consecutive subsequence that is
#'              *mostly* regularly sampled, allowing for a defined number of deviations.
#'              It works by finding all perfectly regular "seed" sequences and then
#'              extending them outwards as long as the sampling interval is within a
#'              given deviation.
#' @param timestamps A numeric vector representing the time of each sample, sorted ascending.
#' @param allowed_deviation The maximum allowed deviation from the common interval
#'                          for a sample to be included in the sequence. For a 7-day
#'                          interval, an `allowed_deviation` of 1 would include
#'                          intervals of 6, 7, and 8.
#' @param min_seed_length The minimum number of perfectly consecutive samples required
#'                        to form a "seed" sequence for extension. Must be at least 2.
#' @param tolerance A small numeric value for floating-point comparisons.
#' @return A list containing details of the best flexible sequence found:
#'         `start_index`, `end_index`, `sequence_length`, `common_interval` (of the seed),
#'         and `mismatch_count` (number of deviated intervals included).
find_longest_flexible_sequence <- function(timestamps, allowed_deviation = 1, min_seed_length = 3, tolerance = 1e-10) {
  
  # --- Step 1: Find all perfectly regular "seed" sequences ---
  if (length(timestamps) < min_seed_length) {
    message("Not enough data points to find a seed sequence.")
    return(NULL)
  }
  intervals <- diff(timestamps)
  if (length(intervals) < 1) return(NULL)
  
  interval_groups <- cumsum(c(1, abs(diff(intervals)) > tolerance))
  runs <- rle(interval_groups)
  
  # Identify all runs that meet the minimum seed length criteria
  seed_indices <- which(runs$lengths >= (min_seed_length - 1))
  
  if (length(seed_indices) == 0) {
    message(paste("No perfect seed sequences of at least length", min_seed_length, "found."))
    return(NULL)
  }
  
  # --- Step 2: Extend each seed and find the best resulting sequence ---
  best_sequence <- list(sequence_length = 0)
  
  for (i in seed_indices) {
    # Determine properties of the current seed
    run_length <- runs$lengths[i]
    run_end_in_intervals <- sum(runs$lengths[1:i])
    run_start_in_intervals <- run_end_in_intervals - run_length + 1
    
    # These are the indices in the original timestamps vector
    seed_start <- run_start_in_intervals
    seed_end <- run_end_in_intervals + 1
    common_interval <- intervals[seed_start]
    
    # --- Extend Forward ---
    extended_end <- seed_end
    mismatch_count <- 0
    
    # Start checking from the sample right after the seed ends
    if (extended_end < length(timestamps)) {
      for (j in (extended_end + 1):length(timestamps)) {
        current_interval <- timestamps[j] - timestamps[j-1]
        
        if (abs(current_interval - common_interval) <= tolerance) {
          # Perfect match, continue
          extended_end <- j
        } else if (abs(current_interval - common_interval) <= allowed_deviation) {
          # Mismatch, but within allowed deviation
          extended_end <- j
          mismatch_count <- mismatch_count + 1
        } else {
          # Deviation is too large, stop extending
          break
        }
      }
    }
    
    # --- Extend Backward ---
    extended_start <- seed_start
    
    # Start checking from the sample right before the seed begins
    if (extended_start > 1) {
      for (j in (extended_start - 1):1) {
        current_interval <- timestamps[j+1] - timestamps[j]
        
        if (abs(current_interval - common_interval) <= tolerance) {
          # Perfect match
          extended_start <- j
        } else if (abs(current_interval - common_interval) <= allowed_deviation) {
          # Mismatch
          extended_start <- j
          mismatch_count <- mismatch_count + 1
        } else {
          # Stop extending
          break
        }
      }
    }
    
    # --- Compare with the best sequence found so far ---
    current_length <- extended_end - extended_start + 1
    if (current_length > best_sequence$sequence_length) {
      best_sequence <- list(
        start_index = extended_start,
        end_index = extended_end,
        sequence_length = current_length,
        common_interval = common_interval,
        mismatch_count = mismatch_count
      )
    }
  }
  
  return(best_sequence)
}


# --- --- --- --- --- ---
#      3. Main Script Logic (Your Custom Code)
# --- --- --- --- --- ---

# Load the initial dataset
d_initial <- readRDS(here("data", "d_initial_Simon_subset.rds"))

# Filter samples from "Randers"
d_filtered <- amp_filter_samples(d_initial, SampleSite == "Randers")

# Create a numeric time series, normalized to start at 0 days.
time_series <- as.numeric(d_filtered$metadata$SampleDate) - as.numeric(d_filtered$metadata$SampleDate[1])

# Ensure the data is sorted
time_series <- sort(time_series)


# --- --- --- --- --- ---
#      4. Run Analysis & Print Results
# --- --- --- --- --- ---

# Find the longest flexible sequence in your time series data.
# We are looking for a sequence with a dominant interval, but allowing for
# intervals that are off by up to 1 day.
consec_seq <- find_longest_flexible_sequence(
  time_series, 
  allowed_deviation = 1, 
  min_seed_length = 2
)

# Print the results
if (length(consec_seq) > 1) {
  cat("--- Longest Flexible Sequence Found in 'Randers' data ---\n")
  cat("Start Index:        ", consec_seq$start_index, "\n")
  cat("End Index:          ", consec_seq$end_index, "\n")
  cat("Total Length:       ", consec_seq$sequence_length, "samples\n")
  cat("Dominant Interval:  ", round(consec_seq$common_interval, 4), "days\n")
  cat("Imperfect Intervals:", consec_seq$mismatch_count, "(off by max", 1, "day)\n\n")
  cat("Span of time series:", time_series[consec_seq$end_index] - time_series[consec_seq$start_index], "days, or", 
      (time_series[consec_seq$end_index] - time_series[consec_seq$start_index])/7, "weeks\n\n")
  
  
  cat("\nIntervals in this sequence (in days):\n")
  print(diff(time_series[consec_seq$start_index:consec_seq$end_index]))
  

  
} else {
  cat("Could not find a suitable flexible sequence in the 'Randers' data.\n")
}

# --- Original Strict Function (kept for reference) ---
# find_longest_regular_sequence <- function(timestamps, tolerance = 1e-10) { ... }

