import os
import pandas as pd
from interpolation_methods import knn_interpolation

# --- 1. CONFIGURE YOUR SCRIPT ---
# TODO: Update these paths to match your folder locations.
# Use raw strings (r"...") to avoid issues with backslashes.
input_folder = r"C:\Users\Simon\Documents\pp_dy\data\abundance_tables"
output_folder = r"C:\Users\Simon\Documents\pp_dy\data\interpolated_abundance_tables"

# The desired interval for your new, regularly sampled data (e.g., 7.0 for every 7 days)
new_time_interval = 7.0

# The number of nearest neighbors to use for the KNN interpolation
K = 5
# --------------------------------

# --- Script execution starts here ---

# Create the output folder if it doesn't exist
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    print(f"Created output folder: {output_folder}")

# Check if the input folder exists
if not os.path.exists(input_folder):
    print(f"ERROR: Input folder not found at {input_folder}")
    print("Please update the 'input_folder' variable in this script.")
else:
    # Get a list of all the .tsv files in the input folder
    files_to_process = [f for f in os.listdir(input_folder) if f.endswith('.tsv')]

    if not files_to_process:
        print(f"No .tsv files found in {input_folder}")
    else:
        # Process each file
        for filename in files_to_process:
            print(f"\nProcessing {filename}...")
            file_path = os.path.join(input_folder, filename)

            # Load the irregularly sampled data from the .tsv file
            irregular_data = pd.read_csv(file_path, sep='\t', index_col=0, low_memory=False)

            # --- Data Cleaning: Keep only columns with valid numeric names ---
            valid_numeric_columns = []
            for col in irregular_data.columns:
                try:
                    float(col)
                    valid_numeric_columns.append(col)
                except ValueError:
                    print(f"  - Ignoring non-numeric column: '{col}'")
            
            irregular_data = irregular_data[valid_numeric_columns]
            # -----------------------------------------------------------------
            
            # Convert column names to float for calculations
            time_points = [float(col) for col in irregular_data.columns]

            if not time_points:
                print(f"  - SKIPPING: No valid time point columns found in {filename}.")
                continue

            # Determine the new, regularly spaced time points
            start_time = min(time_points)
            end_time = max(time_points)
            
            new_time_points = []
            current_time = start_time
            while current_time <= end_time:
                new_time_points.append(current_time)
                current_time += new_time_interval

            # Create a new DataFrame to store the interpolated data
            interpolated_data = pd.DataFrame(index=irregular_data.index)

            # For each new time point, perform KNN interpolation
            for t in new_time_points:
                print(f"  Interpolating for time point: {t:.2f}")
                interpolated_data[t] = knn_interpolation(irregular_data, t, K)
                print

            # Construct the path for the output file
            output_filepath = os.path.join(output_folder, filename)

            # Save the new, regularly sampled DataFrame to a .tsv file
            interpolated_data.to_csv(output_filepath, sep='\t')
            print(f"Finished processing {filename}. Saved to {output_filepath}")

        print("\nAll files have been interpolated!")

