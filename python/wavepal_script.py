# Import necessary libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from wavepal import Wavepal
import os

# --- Configuration ---
# TODO: UPDATE THESE PATHS TO MATCH YOUR FILE LOCATIONS
TIME_SERIES_FILE = '/mnt/c/Users/Simon/Documents/pp_dy/data/abundance_tables_rclr/Randers_abundance_table.tsv'
PAIRINGS_FILE = '/mnt/c/Users/Simon/Documents/pp_dy/output/wavepal_input/Randers_rclr_abund_norm_off_wavepal.txt'
OUTPUT_DIR = '/mnt/c/Users/Simon/Documents/pp_dy/output/wavepal_plots/' # A folder to save the output plots

# --- Main Script ---

def main():
    """
    Main function to load data, iterate through pairs,
    and generate paired time-frequency analysis plots.
    """
    print("Loading data...")
    try:
        ts_df = pd.read_csv(TIME_SERIES_FILE, sep='\t', index_col=0, header=0)
        pairings_df = pd.read_csv(PAIRINGS_FILE, sep='\t')
    except IOError as e:
        print("ERROR: Could not find a file. Please check your file paths. Details: {0}".format(e))
        return
    except Exception as e:
        print("An error occurred during data loading: {0}".format(e))
        return

    numeric_cols = pd.to_numeric(ts_df.columns, errors='coerce')
    valid_cols = ts_df.columns[~np.isnan(numeric_cols)]
    ts_df = ts_df[valid_cols]
    time_points = ts_df.columns.astype(float).values

    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    print("Output plots will be saved to '{0}'".format(OUTPUT_DIR))
    
    theta = np.linspace(time_points[0], time_points[-1], 2000)
    percentile = np.array([95.0])
    ANALYSIS_PARAMS = {
        'theta': theta,
        'w0': 5.5,
        'permin': 10.,
        'percentile': percentile
    }

    for index, row in pairings_df.iterrows():
        var_x = str(row['X'])
        var_y = str(row['Y'])

        print("\nProcessing pair: {0} and {1}...".format(var_x, var_y))

        if var_x not in ts_df.index or var_y not in ts_df.index:
            print("  WARNING: Could not find one or both variables for this pair. Skipping.")
            print("  (Looking for '{0}' and '{1}' in the time series index)".format(var_x, var_y))
            continue

        value_x = ts_df.loc[var_x].values.astype(float)
        value_y = ts_df.loc[var_y].values.astype(float)
        
        print("  Analyzing {0}...".format(var_x))
        wp_x = Wavepal(time_points, value_x) 
        wp_x.check_data()
        wp_x.choose_trend_degree(pol_degree=-1)
        wp_x.trend_vectors()
        wp_x.timefreq_analysis(**ANALYSIS_PARAMS)

        print("  Analyzing {0}...".format(var_y))
        wp_y = Wavepal(time_points, value_y)
        wp_y.check_data()
        wp_y.choose_trend_degree(pol_degree=-1)
        wp_y.trend_vectors()
        wp_y.timefreq_analysis(**ANALYSIS_PARAMS)
        
        if not wp_x.run_timefreq_analysis or not wp_y.run_timefreq_analysis:
            print("  WARNING: Analysis failed for one or both variables. Skipping plot.")
            continue

        print("  Generating plot...")
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 10), sharex=True)

        # --- MODIFIED: Changed y-axis to a logarithmic scale of periods in days ---
        # For plot X
        vmin_x = np.min(wp_x.scalogram)
        vmax_x = np.max(wp_x.scalogram)
        levels_x = np.linspace(vmin_x, vmax_x, 100)
        # 1. Use the raw period values for the y-axis
        axes[0].contourf(wp_x.theta, wp_x.period_cwt, np.transpose(wp_x.scalogram), levels=levels_x, cmap='viridis')
        # 2. Set the scale to logarithmic
        axes[0].set_yscale('log')
        axes[0].set_title("Wavelet Scalogram for: {0}".format(var_x), fontsize=14)
        # 3. Update the axis label
        axes[0].set_ylabel("Period (days)", fontsize=12)
        axes[0].grid(True, linestyle=':', alpha=0.6)

        # For plot Y
        vmin_y = np.min(wp_y.scalogram)
        vmax_y = np.max(wp_y.scalogram)
        levels_y = np.linspace(vmin_y, vmax_y, 100)
        # 1. Use the raw period values for the y-axis
        axes[1].contourf(wp_y.theta, wp_y.period_cwt, np.transpose(wp_y.scalogram), levels=levels_y, cmap='viridis')
        # 2. Set the scale to logarithmic
        axes[1].set_yscale('log')
        axes[1].set_title("Wavelet Scalogram for: {0}".format(var_y), fontsize=14)
        # 3. Update the axis label
        axes[1].set_ylabel("Period (days)", fontsize=12)
        axes[1].set_xlabel("Time (days)", fontsize=12)
        axes[1].grid(True, linestyle=':', alpha=0.6)
        # --- END OF MODIFICATION ---

        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
        fig.suptitle("Time-Frequency Analysis: {0} vs. {1}".format(var_x, var_y), fontsize=16, fontweight='bold')
        
        figname = "{0}_vs_{1}.png".format(var_x, var_y)
        fig.savefig(os.path.join(OUTPUT_DIR, figname), dpi=300)
        plt.close(fig)

    print("\nAll pairs processed successfully!")


if __name__ == '__main__':
    main()