import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

TAXONOMY_COLUMNS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'NCBI_taxonomy_ID']


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


def load_species_table(input_file):
    """Read a taxonomy x samples CSV and collapse it to a species x samples table."""
    df = pd.read_csv(input_file)
    sample_columns = [col for col in df.columns if col not in TAXONOMY_COLUMNS]
    return df.groupby(['species', 'NCBI_taxonomy_ID'])[sample_columns].sum().fillna(0)


def compute_taxa_metrics(df):
    return pd.DataFrame({"percent present": 100 * (df > 0).sum(axis=1) / df.shape[1],
                         "mean reads": df.mean(axis=1),
                         "total reads": df.sum(axis=1)})


def rate_threshold_index(y):
    """Index where the slope of a sorted curve overtakes its average rate of change.

    This is the same criterion as the reads-derivative method (reads lost per species
    removed > the uniform rate), applied to the discrete data instead of a fitted curve.
    Single noisy steps can exceed the average rate early, so the slope is compared on
    balance: the threshold is where the running sum of (average rate - slope) peaks,
    i.e. the last point before which the curve has risen more slowly than average.
    Geometrically this is the knee, the point farthest below the chord joining the ends.

    Returns the index, the discrete slope and the average rate.
    """
    y = np.asarray(y, dtype=float)
    slope = np.diff(y)
    average_rate = (y[-1] - y[0]) / (len(y) - 1)
    rate_deficit = np.concatenate(([0], np.cumsum(average_rate - slope)))
    return int(np.argmax(rate_deficit)), slope, average_rate


def combined_score_threshold(taxa_metrics, save_plots=False, output_dir="."):
    """Keep taxa above the knee of a score combining prevalence and abundance.

    Each taxon is scored as prevalence (fraction of samples) x log10(mean reads + 1).
    The log transform damps the few extremely abundant taxa, and the product means a
    taxon only scores high when it is both prevalent and abundant. Raw magnitudes are
    used rather than ranks, which would give a near-linear curve with no real knee.
    """
    score = taxa_metrics['percent present'] / 100 * np.log10(taxa_metrics['mean reads'] + 1)

    sorted_score = np.sort(score.values)
    knee_idx, score_slope, average_rate = rate_threshold_index(sorted_score)
    score_threshold = sorted_score[knee_idx]
    keep = score > score_threshold

    total_taxa = len(taxa_metrics)
    total_reads = taxa_metrics['mean reads'].sum()
    _stats = {
        'method': 'combined',
        'score_threshold': score_threshold,
        'species_removed': int((~keep).sum()),
        'species_pct_removed': 100 * (~keep).sum() / total_taxa,
        'reads_pct_removed': 100 * taxa_metrics.loc[~keep, 'mean reads'].sum() / total_reads,
        'min_reads_retained': taxa_metrics.loc[keep, 'mean reads'].min(),
        'min_prevalence_retained': taxa_metrics.loc[keep, 'percent present'].min(),
        'total_taxa': total_taxa
    }

    if save_plots:
        # Plot 1: Sorted combined score with knee
        plt.figure(figsize=(10, 6))
        plt.plot(np.arange(total_taxa), sorted_score, label='Combined score')
        plt.axvline(knee_idx, color='black', linestyle='--',
                    label=f'Knee: {_stats["species_removed"]} species removed')
        plt.scatter([knee_idx], [score_threshold], color='black', zorder=5)
        plt.title('Combined Prevalence-Abundance Score')
        plt.xlabel('Number of Lowest-Scoring Species Removed')
        plt.ylabel('Prevalence x log10(mean reads + 1)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "combined_score_knee.png"))
        plt.close()

        # Plot 2: Derivative of the sorted score against its average rate
        plt.figure(figsize=(10, 6))
        plt.plot(np.arange(1, total_taxa), score_slope, linewidth=0.8, label='Score slope (per species)')
        plt.axhline(average_rate, color='tab:red', linestyle=':', label=f'Average rate: {average_rate:.2e}')
        plt.axvline(knee_idx, color='black', linestyle='--',
                    label=f'Threshold: {_stats["species_removed"]} species removed')
        plt.yscale('log')
        plt.title('Derivative of Combined Score as a function of Taxa')
        plt.xlabel('Number of Lowest-Scoring Species Removed')
        plt.ylabel('Rate of Change (log scale)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "combined_score_derivative.png"))
        plt.close()

        # Plot 3: Prevalence vs abundance, coloured by retained/removed
        plt.figure(figsize=(10, 6))
        plt.scatter(taxa_metrics.loc[~keep, 'mean reads'], taxa_metrics.loc[~keep, 'percent present'],
                    s=4, alpha=0.4, color='grey', label=f'Removed ({(~keep).sum()})')
        plt.scatter(taxa_metrics.loc[keep, 'mean reads'], taxa_metrics.loc[keep, 'percent present'],
                    s=4, alpha=0.6, color='tab:blue', label=f'Retained ({keep.sum()})')
        plt.xscale('log')
        plt.title('Species Retained by Combined Score')
        plt.xlabel('Mean Reads (log scale)')
        plt.ylabel('Percent of Samples Present (%)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "combined_score_retained.png"))
        plt.close()

    return keep, _stats


def reads_derivative_threshold(taxa_metrics, save_plots=False, output_dir="."):
    """Keep taxa above the point where the fitted reads-loss rate exceeds the species-loss rate."""
    total_taxa = len(taxa_metrics)

    # Taxa curve
    # Sort by mean reads descending
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
    keep = taxa_metrics['mean reads'] >= min_reads_retained

    # Prepare stats dictionary
    _stats = {
        'method': 'reads',
        'species_removed': int(exceed_idx),
        'species_pct_removed': species_pct_removed,
        'reads_pct_removed': reads_pct_removed,
        'min_reads_retained': min_reads_retained,
        'min_prevalence_retained': taxa_metrics.loc[keep, 'percent present'].min(),
        'species_model_params': popt_log,
        'reads_model_params': popt_sig,
        'total_taxa': total_taxa
    }

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
        plt.xlabel('Number of Lowest-Abundance Species Removed')
        plt.ylabel('Cumulative Percentage (%)')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "removal_curves.png"))
        plt.close()

        # Plot 2: Derivative
        plt.figure(figsize=(10, 6))
        plt.plot(x_reads, reads_derivative, label='Reads Derivative')
        plt.axvline(exceed_idx, color='black', linestyle='--', label=f'Threshold: {exceed_idx} species')
        plt.scatter([exceed_idx], [reads_derivative[exceed_idx]], color='black', zorder=5)

        # Add derivative equation
        deriv_eq = (rf'$\frac{{dy}}{{dx}} = '
                    rf'\frac{{{L_sig:.1f} \cdot {k_sig:.3f} \cdot e^{{-{k_sig:.3f}(x - {x0_sig:.1f})}}}}'
                    rf'{{(1 + e^{{-{k_sig:.3f}(x - {x0_sig:.1f})}})^2}}$')
        plt.text(0.02, 0.98, deriv_eq, transform=plt.gca().transAxes,
                 fontsize=14, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.5))

        plt.title('Derivative of Read Removal Curve as a function of Taxa')
        plt.xlabel('Number of Lowest-Abundance Species Removed')
        plt.ylabel('Rate of Change (%)')
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, "removal_derivatives.png"))
        plt.close()

    return keep, _stats


METHODS = {
    'combined': combined_score_threshold,
    'reads': reads_derivative_threshold,
}


def analyze_normalized_reads(input_file, method='combined', save_plots=False, save_filtered=False, output_dir="."):
    df = load_species_table(input_file)
    taxa_metrics = compute_taxa_metrics(df)

    if save_filtered or save_plots:
        os.makedirs(output_dir, exist_ok=True)

    keep, _stats = METHODS[method](taxa_metrics, save_plots=save_plots, output_dir=output_dir)

    # Retained taxa with their metrics
    filtered_taxa = taxa_metrics[keep].reset_index().rename(columns={'species': 'Taxon'})

    if save_filtered:
        filtered_taxa.to_csv(os.path.join(output_dir, "filtered_taxa.csv"), index=True)

    return _stats, filtered_taxa


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determine filtering threshold based on species prevalence and read abundance.")
    parser.add_argument("input_file", help="Path to normalized read counts CSV file (species x samples)")
    parser.add_argument("--method", choices=METHODS.keys(), default="combined",
                        help="combined: knee of a prevalence x log abundance score (default); "
                             "reads: derivative of the fitted cumulative reads-loss curve")
    parser.add_argument("--save-plots", action="store_true", help="Save output plots (default: no plots)")
    parser.add_argument("--save-filtered", action="store_true", help="Save filtered taxa to CSV file (default: no)")
    parser.add_argument("--output-dir", default=".",
                        help="Directory for saved plots and CSV (default: current directory)")
    args = parser.parse_args()

    stats, filtered_taxa = analyze_normalized_reads(args.input_file,
                                                    method=args.method,
                                                    save_plots=args.save_plots,
                                                    save_filtered=args.save_filtered,
                                                    output_dir=args.output_dir)
    print(f"Method: {stats['method']}")
    if stats['method'] == 'combined':
        print(f"Combined score threshold: {stats['score_threshold']:.4f}")
    print(f"Species removed: {stats['species_removed']} out of {stats['total_taxa']} "
          f"({stats['species_pct_removed']:.2f}%)")
    print(f"Reads removed: {stats['reads_pct_removed']:.2f}%")
    print(f"Minimum mean reads among retained species: {stats['min_reads_retained']:.1f}")
    print(f"Minimum prevalence among retained species: {stats['min_prevalence_retained']:.2f}%")
