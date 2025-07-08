import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def logistic(x, L, k, x0):
    return L / (1 + np.exp(-k * (x - x0)))

def logistic_derivative(x, L, k, x0):
    exp_term = np.exp(-k * (x - x0))
    return (L * k * exp_term) / ((1 + exp_term) ** 2)

def generalized_sigmoid(x, L, k, x0, c):
    return L / (1 + np.exp(-k * (x - x0))) + c

def generalized_sigmoid_derivative(x, L, k, x0, c):
    exp_term = np.exp(-k * (x - x0))
    return (L * k * exp_term) / ((1 + exp_term) ** 2)

def analyze_normalized_reads(input_file, group_by='species', save_plots=False):
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

    # Sort by prevalence (ascending)
    sorted_metrics = taxa_metrics.sort_values(by='percent present').reset_index(drop=True)
    total_species = len(sorted_metrics)

    # Species curve
    x_species = np.arange(1, total_species + 1)
    y_species = np.arange(1, total_species + 1) / total_species * 100
    popt_log, _ = curve_fit(logistic, x_species, y_species, p0=[100, 0.01, total_species/2], maxfev=5000)
    species_fit = logistic(x_species, *popt_log)
    species_derivative = logistic_derivative(x_species, *popt_log)

    # Reads curve
    x_reads = np.arange(0, total_species + 1)
    y_reads_raw = np.concatenate(([0], sorted_metrics['total reads'].cumsum().values))
    y_reads_pct = y_reads_raw / y_reads_raw[-1] * 100
    popt_sig, _ = curve_fit(generalized_sigmoid, x_reads, y_reads_pct, p0=[100, -0.01, total_species/2, 0], maxfev=5000)
    reads_fit = generalized_sigmoid(x_reads, *popt_sig)
    reads_derivative = generalized_sigmoid_derivative(x_reads, *popt_sig)

    # Threshold: when reads loss > species loss per unit
    species_loss_rate = 100 / total_species
    exceed_idx = np.where(reads_derivative > species_loss_rate)[0][0]
    species_pct_removed = (exceed_idx / total_species) * 100
    reads_pct_removed = y_reads_pct[exceed_idx]
    min_reads_retained = sorted_metrics.iloc[exceed_idx:]['total reads'].min()

    _stats = {
        'species_removed': exceed_idx,
        'species_pct_removed': species_pct_removed,
        'reads_pct_removed': reads_pct_removed,
        'min_reads_retained': min_reads_retained,
        'species_model_params': popt_log,
        'reads_model_params': popt_sig
    }

    if save_plots:
        # Plot 1: Species and Reads fits
        plt.figure(figsize=(10, 6))
        plt.plot(x_species, y_species, 'o', markersize=3, alpha=0.6, label='Species Removed (%)')
        plt.plot(x_species, species_fit, label='Fitted Species Logistic')
        plt.plot(x_reads, y_reads_pct, 'o', markersize=3, alpha=0.6, label='Reads Removed (%)')
        plt.plot(x_reads, reads_fit, label='Fitted Reads Sigmoid')
        plt.axvline(exceed_idx, color='black', linestyle='--', label=f'Threshold: {exceed_idx} species')
        plt.legend()
        plt.title('Species and Reads Removal Curves')
        plt.xlabel('Number of Least Prevalent Species Removed')
        plt.ylabel('Cumulative Percentage (%)')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("removal_curves.png")
        plt.close()

        # Plot 2: Derivatives
        plt.figure(figsize=(10, 6))
        plt.plot(x_species, species_derivative, label='Species Derivative')
        plt.plot(x_reads, reads_derivative, label='Reads Derivative')
        plt.axhline(species_loss_rate, color='blue', linestyle='--', label='Species Loss Rate')
        plt.axvline(exceed_idx, color='black', linestyle='--', label=f'Threshold: {exceed_idx} species')
        plt.scatter([exceed_idx], [reads_derivative[exceed_idx]], color='black', zorder=5)
        plt.title('Derivatives of Removal Curves')
        plt.xlabel('Number of Least Prevalent Species Removed')
        plt.ylabel('Rate of Change (%)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig("removal_derivatives.png")
        plt.close()

    return exceed_idx, _stats

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Determine filtering threshold based on species prevalence and read abundance.")
    parser.add_argument("input_file", help="Path to normalized read counts CSV file (species x samples)")
    parser.add_argument("--save-plots", action="store_true", help="Save output plots (default: no plots)")
    args = parser.parse_args()

    threshold, stats = analyze_normalized_reads(args.input_file, save_plots=args.save_plots)
    print(f"Recommended filtering threshold: {threshold} species removed")
    print(f"Species removed: {stats['species_pct_removed']:.2f}%, Reads removed: {stats['reads_pct_removed']:.2f}%")
    print(f"Minimum reads among retained species: {stats['min_reads_retained']:.1f}")
