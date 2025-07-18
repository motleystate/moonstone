import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Define the bounded logistic model (taxa model)
def bounded_logistic(x, L, a, b):
    return L / (1 + a * x**b)

# Define generalized sigmoid function (reads model)
def generalized_sigmoid(x, L, k, x0, c):
    return L / (1 + np.exp(-k * (x - x0))) + c

# Define derivative of the generalized sigmoid function
def generalized_sigmoid_derivative(x, L, k, x0, c):
    exp_term = np.exp(-k * (x - x0))
    return (L * k * exp_term) / ((1 + exp_term) ** 2)

def analyze_normalized_reads(input_file, group_by='species', save_plots=False, save_filtered=False):
    df = pd.read_csv(input_file)
    # Melt the dataframe
    melted = df.melt(
        id_vars=['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'NCBI_taxonomy_ID'],
        value_vars=[col for col in df.columns if col not in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'NCBI_taxonomy_ID']],
        var_name='Sample',
        value_name='Abundance'
    )
    # Group by species and taxonomy ID along with Sample
    grouped = melted.groupby(['species', 'NCBI_taxonomy_ID', 'Sample'], as_index=False)['Abundance'].sum()

    # Pivot with a MultiIndex on rows: species + NCBI_taxonomy_ID
    df = grouped.pivot(index=['species', 'NCBI_taxonomy_ID'], columns='Sample', values='Abundance').fillna(0)

    taxa_metrics = pd.DataFrame({"percent present": 100 * (1 - (df == 0).astype(int).sum(axis=1) / df.shape[1]),
                                "mean reads": df.mean(axis=1),
                                "total reads": df.sum(axis=1)})
    
    total_taxa = len(taxa_metrics)

    # Taxa curve
    # Sort by percent present descending
    sorted_metrics = taxa_metrics.sort_values(by='mean reads', ascending=False).reset_index(drop=True)
    x_taxa = np.arange(1, len(sorted_metrics) + 1)
    y_taxa = sorted_metrics['percent present'].values
    
    # Set dynamic bounds for L based on data max, fit and predict
    max_present = y_taxa.max()
    bounds_dynamic = ([max_present - 10, 0, -np.inf], [max_present + 10, np.inf, np.inf])
    popt_log, _ = curve_fit(bounded_logistic, x_taxa, y_taxa, bounds=bounds_dynamic, maxfev=5000)
    taxa_fit = bounded_logistic(x_taxa, *popt_log)

    # Calculate Taxa Curve R²
    ss_res = np.sum((y_taxa - taxa_fit) ** 2)
    ss_tot = np.sum((y_taxa - np.mean(y_taxa)) ** 2)
    taxa_r_squared = 1 - (ss_res / ss_tot)
    
    # Reads curve
    sorted_metrics = taxa_metrics.sort_values(by='mean reads').reset_index(drop=True)
    x_reads = np.arange(0, total_taxa + 1)
    y_reads_raw = np.concatenate(([0], sorted_metrics['mean reads'].cumsum().values))
    y_reads_pct = y_reads_raw / y_reads_raw[-1] * 100
    popt_sig, _ = curve_fit(generalized_sigmoid, x_reads, y_reads_pct, p0=[100, -0.01, total_taxa/2, 0], maxfev=5000)
    reads_fit = generalized_sigmoid(x_reads, *popt_sig)
    
    # Calculate Read Curve R²    
    ss_res_final = np.sum((y_reads_pct - reads_fit) ** 2)
    ss_tot_final = np.sum((y_reads_pct - np.mean(y_reads_pct)) ** 2)
    reads_r_squared = 1 - (ss_res_final / ss_tot_final)
    
    # Calculate Derivative of read sigmoid function
    reads_derivative = generalized_sigmoid_derivative(x_reads, *popt_sig)

    # Threshold: when reads loss > species loss per unit
    species_loss_rate = 100 / total_taxa
    exceed_idx = np.where(reads_derivative > species_loss_rate)[0][0]
    species_pct_removed = (exceed_idx / total_taxa) * 100
    reads_pct_removed = y_reads_pct[exceed_idx]
    min_reads_retained = sorted_metrics.iloc[exceed_idx:]['mean reads'].min()

    # Prepare stats dictionary
    _stats = {
        'species_removed': exceed_idx,
        'species_pct_removed': species_pct_removed,
        'reads_pct_removed': reads_pct_removed,
        'min_reads_retained': min_reads_retained,
        'species_model_params': popt_log,
        'reads_model_params': popt_sig,
        'total_taxa': total_taxa
    }
    # Make a filtered version of the `taxa_metrics` DataFrame
    filtered_taxa = pd.DataFrame({"Taxon": df.index.get_level_values('species'),
                                "NCBI_taxonomy_ID": df.index.get_level_values('NCBI_taxonomy_ID'),
                                "percent present": 100 * (1 - (df == 0).astype(int).sum(axis=1) / df.shape[1]),
                                "mean reads": df.mean(axis=1),
                                "total reads": df.sum(axis=1)})
    filtered_taxa = filtered_taxa.query('`mean reads` >= @_stats["min_reads_retained"]').reset_index(drop=True)

    if save_filtered:
        filtered_taxa.to_csv("filtered_taxa.csv", index=True)


    if save_plots:
        # Plot 1: Species and Reads fits
        plt.figure(figsize=(10, 6))
        plt.plot(x_taxa, y_taxa, 'o', markersize=3, alpha=0.6, label='Species (%)')
        plt.plot(x_taxa, taxa_fit, label=f'Fitted Species Logistic (R² = {taxa_r_squared:.4f})')
        plt.plot(x_reads, y_reads_pct, 'o', markersize=3, alpha=0.6, label='Reads (%)')
        plt.plot(x_reads, reads_fit, label=f'Fitted Reads Sigmoid (R² = {reads_r_squared:.4f})')
        plt.legend()

         # Add equations as text
        L_log, a_log, b_log = popt_log
        taxa_eq = rf'$y = \frac{{{L_log:.1f}}}{{1 + {a_log:.2e} \cdot x^{{{b_log:.2f}}}}}$'
        
        L_sig, k_sig, x0_sig, c_sig = popt_sig
        reads_eq = rf'$y = \frac{{{L_sig:.1f}}}{{1 + e^{{-{k_sig:.3f}(x - {x0_sig:.1f})}}}} + {c_sig:.2f}$'
        
        # Position equations in upper right corner
        plt.text(0.98, 0.98, 'Species: ' + taxa_eq, transform=plt.gca().transAxes, 
                fontsize=12, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        plt.text(0.98, 0.88, 'Reads: ' + reads_eq, transform=plt.gca().transAxes, 
                fontsize=12, verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))

        plt.title('Species and Reads Curves')
        plt.xlabel('Number of Least Prevalent Species Removed')
        plt.ylabel('Cumulative Percentage (%)')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("removal_curves.png")
        plt.close()

        # Plot 2: Derivative
        plt.figure(figsize=(10, 6))
        plt.plot(x_reads, reads_derivative, label='Reads Derivative')
        plt.axvline(exceed_idx, color='black', linestyle='--', label=f'Threshold: {exceed_idx} species')
        plt.scatter([exceed_idx], [reads_derivative[exceed_idx]], color='black', zorder=5)

        # Add derivative equation
        deriv_eq = rf'$\frac{{dy}}{{dx}} = \frac{{{L_sig:.1f} \cdot {k_sig:.3f} \cdot e^{{-{k_sig:.3f}(x - {x0_sig:.1f})}}}}{{(1 + e^{{-{k_sig:.3f}(x - {x0_sig:.1f})}})^2}}$'
        plt.text(0.02, 0.98, deriv_eq, transform=plt.gca().transAxes, 
                fontsize=14, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

        plt.title('Derivative of Read Removal Curve as a function of Taxa')
        plt.xlabel('NNumber of Low-Abundance Species Removed')
        plt.ylabel('Rate of Change (%)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("removal_derivatives.png")
        plt.close()

    return exceed_idx, _stats, filtered_taxa

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine filtering threshold based on species prevalence and read abundance.")
    parser.add_argument("input_file", help="Path to normalized read counts CSV file (species x samples)")
    parser.add_argument("--save-plots", action="store_true", help="Save output plots (default: no plots)")
    parser.add_argument("--save-filtered", action="store_true", help="Save filtered taxa to CSV file (default: no)")
    args = parser.parse_args()

    threshold, stats, filtered_taxa = analyze_normalized_reads(args.input_file, 
                                                               save_plots=args.save_plots,
                                                               save_filtered=args.save_filtered)
    print(f"Recommended filtering threshold: {stats['min_reads_retained']:.1f}")
    print(f"Species removed: {stats['species_removed']} out of {stats['total_taxa']} ({threshold/stats['total_taxa'] * 100:.2f}%)")
    print(f"Reads removed: {stats['reads_pct_removed']:.2f}%")
    print(f"Minimum mean reads among retained species: {stats['min_reads_retained']:.1f}")

