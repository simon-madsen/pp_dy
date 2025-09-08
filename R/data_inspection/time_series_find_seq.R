# --- --- --- --- --- ---
#      1. Load Libraries
# --- --- --- --- --- ---
library(ampvis2)
library(here)


# --- --- --- --- --- ---
#      2. Function Definitions
# --- --- --- --- --- ---

#' @title Find Longest, Most Flexible Consecutive Sequence by Time Span
#' @description This function identifies the consecutive subsequence that spans the
#'              longest amount of time. It is robust to minor deviations, can
#'              "skip" outlier samples, and allows for large "gaps" in the data
#'              based on a set of rules.
#' @param timestamps A numeric vector of timestamps, sorted ascending.
#' @param max_interval The maximum allowed dominant sampling interval for a seed
#'                     sequence to be considered (e.g., 10 days).
#' @param allowed_deviation The maximum allowed deviation from the dominant interval
#'                          for a sample to be included in the sequence.
#' @param min_seed_length The minimum number of perfectly consecutive samples required
#'                        to form a "seed" sequence for extension.
#' @param tolerance A small numeric value for floating-point comparisons.
#' @param gap_allowance_period The number of days a sequence must span to earn one gap.
#' @param gap_extension_threshold The minimum number of days a sequence must continue
#'                                after a gap for the gap to be considered valid.
#' @return A list containing details of the best flexible sequence found, including
#'         indices, lengths, counts for mismatches, skips, and gaps, and the total time span.
find_longest_flexible_sequence <- function(timestamps, max_interval = 10, allowed_deviation = 1, min_seed_length = 3, tolerance = 1e-10, gap_allowance_period = 90, gap_extension_threshold = 21) {
  
  # --- Step 1: Find all perfectly regular "seed" sequences ---
  if (length(timestamps) < min_seed_length) {
    message("Not enough data points to find a seed sequence.")
    return(NULL)
  }
  intervals <- diff(timestamps)
  if (length(intervals) < 1) return(NULL)
  
  interval_groups <- cumsum(c(1, abs(diff(intervals)) > tolerance))
  runs <- rle(interval_groups)
  
  potential_seed_indices <- which(runs$lengths >= (min_seed_length - 1))
  
  if (length(potential_seed_indices) == 0) {
    message(paste("No perfect seed sequences of at least length", min_seed_length, "found."))
    return(NULL)
  }
  
  valid_seed_indices <- c()
  for(i in potential_seed_indices) {
    run_length <- runs$lengths[i]
    run_end_in_intervals <- sum(runs$lengths[1:i])
    run_start_in_intervals <- run_end_in_intervals - run_length + 1
    common_interval <- intervals[run_start_in_intervals]
    if (common_interval <= max_interval) {
      valid_seed_indices <- c(valid_seed_indices, i)
    }
  }
  
  if (length(valid_seed_indices) == 0) {
    message(paste("No seed sequences found with a dominant interval of", max_interval, "days or less."))
    return(NULL)
  }
  
  # --- Step 2: Extend each valid seed and find the best resulting sequence by time span ---
  best_sequence <- list(time_span = 0)
  
  for (i in valid_seed_indices) {
    run_length <- runs$lengths[i]
    run_end_in_intervals <- sum(runs$lengths[1:i])
    run_start_in_intervals <- run_end_in_intervals - run_length + 1
    
    seed_start_idx <- run_start_in_intervals
    seed_end_idx <- run_end_in_intervals + 1
    common_interval <- intervals[seed_start_idx]
    
    # --- Helper function for lookahead/lookbehind validation ---
    validate_extension <- function(start_pos, direction = "forward") {
      if (start_pos < 1 || start_pos > length(timestamps)) return(FALSE)
      
      end_pos <- start_pos
      current_pos <- start_pos
      
      while (TRUE) {
        next_pos <- if (direction == "forward") current_pos + 1 else current_pos - 1
        if (next_pos < 1 || next_pos > length(timestamps)) break
        
        interval <- if (direction == "forward") timestamps[next_pos] - timestamps[current_pos] else timestamps[current_pos] - timestamps[next_pos]
        
        if (abs(interval - common_interval) <= allowed_deviation) {
          end_pos <- next_pos
          current_pos <- next_pos
          span <- if (direction == "forward") timestamps[end_pos] - timestamps[start_pos] else timestamps[start_pos] - timestamps[end_pos]
          if (span >= gap_extension_threshold) return(TRUE)
        } else {
          break
        }
      }
      span <- if (direction == "forward") timestamps[end_pos] - timestamps[start_pos] else timestamps[start_pos] - timestamps[end_pos]
      return(span >= gap_extension_threshold)
    }
    
    # --- Extend Forward ---
    final_end <- seed_end_idx
    forward_skipped <- c(); forward_gaps <- c(); gaps_used_fwd <- 0
    current_pos <- seed_end_idx
    
    while (current_pos < length(timestamps)) {
      next_pos <- current_pos + 1
      if (next_pos > length(timestamps)) break
      
      interval <- timestamps[next_pos] - timestamps[current_pos]
      if (abs(interval - common_interval) <= allowed_deviation) { final_end <- next_pos; current_pos <- next_pos; next }
      
      skip_pos <- next_pos + 1
      if (interval < (common_interval - allowed_deviation) && skip_pos <= length(timestamps)) {
        if (abs((timestamps[skip_pos] - timestamps[current_pos]) - common_interval) <= allowed_deviation) {
          forward_skipped <- c(forward_skipped, next_pos); final_end <- skip_pos; current_pos <- skip_pos; next
        }
      }
      
      span <- timestamps[final_end] - timestamps[seed_start_idx]
      if (gaps_used_fwd < floor(span / gap_allowance_period)) {
        found_gap <- FALSE
        for (j in (next_pos + 1):length(timestamps)) {
          gap_interval <- timestamps[j] - timestamps[current_pos]
          n_mult <- round(gap_interval / common_interval)
          if (n_mult > 1 && n_mult <= 3 && abs(gap_interval - n_mult * common_interval) <= allowed_deviation && validate_extension(j, "forward")) {
            gaps_used_fwd <- gaps_used_fwd + 1; forward_gaps <- c(forward_gaps, (current_pos + 1):(j - 1))
            final_end <- j; current_pos <- j; found_gap <- TRUE; break
          }
        }
        if (found_gap) next
      }
      break
    }
    
    # --- Extend Backward ---
    final_start <- seed_start_idx
    backward_skipped <- c(); backward_gaps <- c(); gaps_used_bwd <- 0
    current_pos <- seed_start_idx
    
    while (current_pos > 1) {
      prev_pos <- current_pos - 1
      if (prev_pos < 1) break
      
      interval <- timestamps[current_pos] - timestamps[prev_pos]
      if (abs(interval - common_interval) <= allowed_deviation) { final_start <- prev_pos; current_pos <- prev_pos; next }
      
      skip_pos <- prev_pos - 1
      if (interval < (common_interval - allowed_deviation) && skip_pos >= 1) {
        if (abs((timestamps[current_pos] - timestamps[skip_pos]) - common_interval) <= allowed_deviation) {
          backward_skipped <- c(backward_skipped, prev_pos); final_start <- skip_pos; current_pos <- skip_pos; next
        }
      }
      
      span <- timestamps[final_end] - timestamps[final_start]
      if ((gaps_used_fwd + gaps_used_bwd) < floor(span / gap_allowance_period)) {
        found_gap <- FALSE
        for (j in (prev_pos - 1):1) {
          gap_interval <- timestamps[current_pos] - timestamps[j]
          n_mult <- round(gap_interval / common_interval)
          if (n_mult > 1 && n_mult <= 3 && abs(gap_interval - n_mult * common_interval) <= allowed_deviation && validate_extension(j, "backward")) {
            gaps_used_bwd <- gaps_used_bwd + 1; backward_gaps <- c(backward_gaps, (j + 1):(current_pos - 1))
            final_start <- j; current_pos <- j; found_gap <- TRUE; break
          }
        }
        if (found_gap) next
      }
      break
    }
    
    # --- Calculate metrics and compare ---
    current_time_span <- timestamps[final_end] - timestamps[final_start]
    if (current_time_span > best_sequence$time_span) {
      total_skipped <- c(backward_skipped, forward_skipped)
      total_gaps <- c(backward_gaps, forward_gaps)
      all_indices <- final_start:final_end
      final_indices <- all_indices[!all_indices %in% c(total_skipped, total_gaps)]
      final_intervals <- diff(timestamps[final_indices])
      
      best_sequence <- list(
        start_index = final_start, end_index = final_end,
        time_span = current_time_span,
        sequence_length = length(final_indices),
        common_interval = common_interval,
        mismatch_count = sum(abs(final_intervals - common_interval) > tolerance),
        skipped_count = length(total_skipped),
        gap_count = gaps_used_fwd + gaps_used_bwd,
        final_indices = final_indices
      )
    }
  }
  
  return(best_sequence)
}


# --- --- --- --- --- ---
#      3. Main Script Logic (Your Custom Code)
# --- --- --- --- --- ---

d_initial <- readRDS(here("data", "d_initial_Simon_subset.rds"))
d_filtered <- amp_filter_samples(d_initial, SampleSite == "Aalborg E")
time_series <- as.numeric(d_filtered$metadata$SampleDate) - as.numeric(d_filtered$metadata$SampleDate[1])
time_series <- sort(time_series)


# --- --- --- --- --- ---
#      4. Run Analysis & Print Results
# --- --- --- --- --- ---

consec_seq <- find_longest_flexible_sequence(
  time_series, 
  max_interval = 10,
  allowed_deviation = 2, 
  min_seed_length = 3,
  gap_allowance_period = 60,
  gap_extension_threshold = 21
)

if (length(consec_seq) > 1) {
  cat("--- Longest Time Span Sequence Found in 'Randers' data ---\n")
  cat("Start Index:        ", consec_seq$start_index, "\n")
  cat("End Index:          ", consec_seq$end_index, "\n")
  cat("Total Time Span:    ", round(consec_seq$time_span, 2), "days\n")
  cat("Samples Used:       ", consec_seq$sequence_length, "\n")
  cat("Dominant Interval:  ", round(consec_seq$common_interval, 4), "days (<= 10 day constraint)\n")
  cat("Imperfect Intervals:", consec_seq$mismatch_count, "(off by max allowed_deviation", ")\n")
  cat("Skipped Samples:    ", consec_seq$skipped_count, "\n")
  cat("Gaps Used:          ", consec_seq$gap_count, "(earned 1 per", 90, "days)\n\n")
  
  #sequence_metadata <- d_filtered$metadata[consec_seq$start_index:consec_seq$end_index, ]
  #cat("Metadata for the full sequence time span:\n")
  #print(head(sequence_metadata))
  
  cat("\nIntervals in this sequence (in days, from non-skipped/gap samples):\n")
  print(diff(time_series[consec_seq$final_indices]))
  
} else {
  cat("Could not find a suitable sequence in the 'Randers' data.\n")
}

