#!/bin/bash
#
# run_elsa.sh
#
# This script runs the eLSA (Extended Local Similarity Analysis) on all
# prepared data tables. It processes files serially, one by one, using the
# lsa_compute command.
#
# It assumes you have an active conda environment where the 'elsa' package
# is installed.
#
# --- Configuration ---

# !!! IMPORTANT: PASTE THE FULL PATH TO lsa_compute HERE !!!
# Find this path by running 'conda activate elsa_env' and then 'which lsa_compute'
LSA_COMPUTE_EXECUTABLE="/home/simon/miniconda3/envs/elsa_env/bin/lsa_compute"


# The script expects to be run from the root of your project directory.
# `.` refers to the current directory.
PROJECT_ROOT="."

# Directory where the R script saved the formatted tables.
INPUT_DIR="$PROJECT_ROOT/output/elsa_tables"

# Directory where the final eLSA results will be saved.
OUTPUT_DIR="$PROJECT_ROOT/output/elsa_results"

# --- eLSA Parameters ---
# -r: Number of replicates per time point. Your format is "tXr1", so you have 1.
REPLICATES=1
# -p: P-value estimation method. "perm" uses permutation.
PVALUE_METHOD="perm"
# -x: Number of permutations for the p-value test.
NUM_PERMUTATIONS=1000
# -d: Maximum time delay to consider. This is a critical parameter.
MAX_DELAY=5
# -n: The normalization method to apply.
NORM_METHOD="percentileZ"

# --- Script Execution ---

# Exit immediately if a command fails
set -e

echo "Starting eLSA computation..."

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"
echo "Results will be saved in: $OUTPUT_DIR"

# Check if input directory is empty
if [ -z "$(ls -A "$INPUT_DIR")" ]; then
   echo "Error: Input directory '$INPUT_DIR' is empty. Did the R script run correctly?"
   exit 1
fi

# Loop through every .tsv file in the input directory
for input_file in "$INPUT_DIR"/*.tsv; do

    # Get the base name of the file (e.g., "Aalborg_W_rclr_abund")
    base_name=$(basename "$input_file" .tsv)

    # Get the number of time points (-s) by counting columns in the header and subtracting the first (#OTU_ID) column.
    num_timepoints=$(head -n 1 "$input_file" | awk -F'\t' '{print NF-1}')

    # Define the output file path
    output_file="$OUTPUT_DIR/${base_name}.lsa"

    echo "----------------------------------------------------"
    echo "Processing: $base_name (Timepoints: $num_timepoints)"

    # Run the lsa_compute command with the correct syntax:
    # lsa_compute <inputFile> <outputFile> [options]
    "$LSA_COMPUTE_EXECUTABLE" "$input_file" "$output_file" \
        -r "$REPLICATES" \
        -p "$PVALUE_METHOD" \
        -d "$MAX_DELAY" \
        -s "$num_timepoints" \
        -x "$NUM_PERMUTATIONS" \
        -n "$NORM_METHOD"

    echo "Finished processing $base_name"
done

echo "----------------------------------------------------"
echo "All analyses are complete."

