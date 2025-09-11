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
OUTPUT_DIR = '/mnt/c/Users/Simon/Documents/pp_dy/output/wavepal_freq_plots/' # A folder to save the output plots

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
    
    #theta = np.linspace(time_points[0], time_points[-1], 2000)
    percentile = np.array([95.0])
    ANALYSIS_PARAMS = {
        #'computes_amplitude': False,
        'D': 750.,
        'n_moments': 10,
        'freqstep': 0.00001,
        'percentile': percentile,
        'weighted_WOSA': True,
        #'freqmin': 1/1000.,
        'freqmax': 1/15.,
        #'WOSA_segments': "all",
        #'freq_min_bound': False
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
        
       # --- MODIFIED: Added carma_params() to enable confidence level calculation ---
        print("  Analyzing {0}...".format(var_x))
        print((time_points[-1]-time_points[0]))
        print(1/(time_points[-1]-time_points[0]))
        wp_x = Wavepal(time_points, value_x) 
        wp_x.check_data()
        wp_x.choose_trend_degree(pol_degree=-1)
        wp_x.trend_vectors()
        wp_x.carma_params(p=1, q=0, signif_level_type='a') # Analytical confidence levels for CAR(1,0) noise
        wp_x.freq_analysis(**ANALYSIS_PARAMS)

        print("  Analyzing {0}...".format(var_y))
        wp_y = Wavepal(time_points, value_y)
        wp_y.check_data()
        wp_y.choose_trend_degree(pol_degree=-1)
        wp_y.trend_vectors()
        wp_y.carma_params(p=1, q=0, signif_level_type='a')
        wp_y.freq_analysis(**ANALYSIS_PARAMS)
        # --- END OF MODIFICATION ---
        
        if not wp_x.run_freq_analysis or not wp_y.run_freq_analysis:
            print("  WARNING: Analysis failed for one or both variables. Skipping plot.")
            continue

        print("  Generating plot...")
        fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 10), sharex=True)

        # --- MODIFIED: Added confidence level plotting ---
        # For plot X
        axes[0].plot(1. / wp_x.freq, wp_x.periodogram, color='black', label='Periodogram')
        # Plot the first (and only) calculated confidence level
        axes[0].plot(1. / wp_x.freq, wp_x.periodogram_cl_anal[:, 0], color='red', linestyle='--', label='95% Confidence Level')
        axes[0].set_title("Periodogram for: {0}".format(var_x), fontsize=14)
        axes[0].set_ylabel("Power", fontsize=12)
        #axes[0].set_xlim([15, 370])
        axes[0].grid(True, linestyle=':', alpha=0.6)
        axes[0].legend()

        # For plot Y
        axes[1].plot(1. / wp_y.freq, wp_y.periodogram, color='black', label='Periodogram')
        axes[1].plot(1. / wp_y.freq, wp_y.periodogram_cl_anal[:, 0], color='red', linestyle='--', label='95% Confidence Level')
        axes[1].set_title("Periodogram for: {0}".format(var_y), fontsize=14)
        axes[1].set_ylabel("Power", fontsize=12)
        axes[1].set_xlabel("Period (days)", fontsize=12)
        axes[1].grid(True, linestyle=':', alpha=0.6)
        axes[1].legend()
        # --- END OF MODIFICATION ---

        plt.tight_layout(rect=[0, 0.03, 1, 0.96])
        fig.suptitle("WOSA Periodogram Analysis: {0} vs. {1}".format(var_x, var_y), fontsize=16, fontweight='bold')
        
        figname = "{0}_vs_{1}_periodogram.png".format(var_x, var_y)
        fig.savefig(os.path.join(OUTPUT_DIR, figname), dpi=300)
        plt.close(fig)
    print("\nAll pairs processed successfully!")


if __name__ == '__main__':
    main()