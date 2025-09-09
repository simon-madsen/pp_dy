import pandas as pd
import numpy as np

# Important: Assumes 'interpolation_methods.py' is in the same directory
from interpolation_methods import knn_interpolation

# --- 1. Generate Sample Data ---
# We create a simple dataset with two microbes and four irregular time points.
# The data is stored in a dictionary and then converted to a pandas DataFrame.
data = {
    '0': [10, 100],   # Time point 0
    '3': [15, 80],    # Time point 3
    '8': [30, 40],    # Time point 8
    '15': [50, 20]    # Time point 15
}
# The index represents the different microbes/features.
index = ['Microbe A', 'Microbe B']
sample_data = pd.DataFrame(data, index=index)

print("--- Original Data ---")
print(sample_data)
print("\n" + "="*30 + "\n")


# --- 2. Set Interpolation Parameters ---
# We want to find the abundance values at a time point where we have no data.
time_to_interpolate = 5.0

# We will use the 2 nearest neighbors to perform the calculation.
K = 2


# --- 3. Run the Interpolation ---
# We call the function from your library file.
# It will return a new column of data for our target time point.
print(f"Interpolating data for time point: {time_to_interpolate} using K={K}...")
interpolated_column = knn_interpolation(sample_data, time_to_interpolate, K)


# --- 4. Display the Results ---
# We add the new, interpolated column to our original DataFrame to see the result.
# Note: The function returns a NumPy array, so we convert it to a pandas Series.
sample_data[str(time_to_interpolate)] = pd.Series(interpolated_column, index=sample_data.index)

# Sort the columns by time point to make the final table easy to read.
final_data = sample_data.reindex(sorted(sample_data.columns, key=float), axis=1)

print("\n--- Data with Interpolated Column ---")
print(final_data)
