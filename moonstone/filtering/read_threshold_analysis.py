import argparse
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

TAXONOMY_COLUMNS = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'NCBI_taxonomy_ID']

# Plot styling shared by all figures
COLORS = {
    'surface': '#fcfcfb',
    'text': '#0b0b0b',
    'text_secondary': '#52514e',
    'grid': '#e6e5e1',
    'primary': '#2a78d6',
    'secondary': '#eb6834',
    'removed': '#b9b7b0',
}
PLOT_STYLE = {
    'figure.facecolor': COLORS['surface'],
    'axes.facecolor': COLORS['surface'],
    'savefig.facecolor': COLORS['surface'],
    'savefig.dpi': 150,
    'font.size': 10,
    'mathtext.fontset': 'cm',
    'text.color': COLORS['text'],
    'axes.edgecolor': COLORS['text_secondary'],
    'axes.labelcolor': COLORS['text'],
    'axes.labelsize': 10.5,
    'axes.titlesize': 14,
    'axes.titleweight': 'bold',
    'axes.titlelocation': 'left',
    'axes.titlecolor': COLORS['text'],
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.grid': True,
    'axes.axisbelow': True,
    'grid.color': COLORS['grid'],
    'grid.linewidth': 0.8,
    'xtick.color': COLORS['text_secondary'],
    'ytick.color': COLORS['text_secondary'],
    'legend.frameon': False,
    'legend.fontsize': 9.5,
    'lines.linewidth': 2,
}


def _num(n):
    """Integer with thousands separators, for use inside mathtext."""
    return f'{n:,}'.replace(',', '{,}')


def _sci(v, digits=2):
    """Number in scientific notation, for use inside mathtext."""
    mantissa, exponent = f'{v:.{digits}e}'.split('e')
    return rf'{mantissa} \times 10^{{{int(exponent)}}}'


def _titles(ax, title, subtitle):
    ax.set_title(title, pad=26)
    ax.text(0, 1.02, subtitle, transform=ax.transAxes, ha='left', va='bottom',
            fontsize=10, color=COLORS['text_secondary'])


def _equation_box(ax, text, x=0.02, y=0.96, ha='left', va='top', fontsize=12):
    ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=fontsize, color=COLORS['text'],
            linespacing=1.6, zorder=6,
            bbox=dict(boxstyle='round,pad=0.6', facecolor='white', edgecolor=COLORS['grid'], linewidth=1))


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

    The criterion is that removing one more taxon changes the curve faster than the
    uniform rate, computed on the discrete data rather than on a fitted curve.
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
        'score_threshold': score_threshold,
        'species_removed': int((~keep).sum()),
        'species_pct_removed': 100 * (~keep).sum() / total_taxa,
        'reads_pct_removed': 100 * taxa_metrics.loc[~keep, 'mean reads'].sum() / total_reads,
        'min_reads_retained': taxa_metrics.loc[keep, 'mean reads'].min(),
        'min_prevalence_retained': taxa_metrics.loc[keep, 'percent present'].min(),
        'total_taxa': total_taxa
    }

    if save_plots:
        plot_combined_score(taxa_metrics, keep, sorted_score, knee_idx, score_slope, average_rate, _stats, output_dir)

    return keep, _stats


def plot_combined_score(taxa_metrics, keep, sorted_score, knee_idx, score_slope, average_rate, _stats, output_dir):
    total_taxa = len(sorted_score)
    # Position k on the x axis is the number of lowest-scoring species removed
    x = np.arange(1, total_taxa + 1)
    k_star = knee_idx + 1
    s_star = sorted_score[knee_idx]
    n_kept, n_removed = int(keep.sum()), int((~keep).sum())
    summary = (f'Keeps {n_kept:,} of {total_taxa:,} taxa ({100 - _stats["species_pct_removed"]:.1f}%), '
               f'holding {100 - _stats["reads_pct_removed"]:.1f}% of reads')
    score_eq = r'$s_i = p_i \, \log_{10}(\bar{r}_i + 1)$'
    score_legend = ('$p_i$: fraction of samples where taxon $i$ is present\n'
                    r'$\bar{r}_i$: mean reads of taxon $i$')

    with plt.rc_context(PLOT_STYLE):
        # Plot 1: Sorted combined score, its chord (average rate) and the threshold
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.axvspan(k_star, total_taxa, color=COLORS['primary'], alpha=0.07, linewidth=0)
        ax.plot([1, total_taxa], [sorted_score[0], sorted_score[-1]], color=COLORS['secondary'],
                linewidth=1.5, linestyle=':', label=r'Chord: average rate $\bar{m}$')
        ax.plot(x, sorted_score, color=COLORS['primary'], label='Sorted score $s_k$')
        ax.axvline(k_star, color=COLORS['text_secondary'], linewidth=1, linestyle='--')
        ax.plot([k_star], [s_star], 'o', markersize=9, color=COLORS['primary'],
                markeredgecolor=COLORS['surface'], markeredgewidth=2, zorder=5)
        ax.annotate(f'$k^* = {_num(k_star)}$\n$s^* = {s_star:.3f}$', xy=(k_star, s_star),
                    xytext=(-16, 40), textcoords='offset points', ha='right', fontsize=12,
                    arrowprops=dict(arrowstyle='-', color=COLORS['text_secondary'], linewidth=0.8))
        ax.text((k_star + total_taxa) / 2, 0.5, f'Retained\n{n_kept:,}', transform=ax.get_xaxis_transform(),
                ha='center', va='center', fontsize=9.5, color=COLORS['text_secondary'])
        _equation_box(ax, score_eq + '\n' + score_legend, fontsize=11)
        _titles(ax, 'Combined Prevalence-Abundance Score', summary)
        ax.set_xlabel('Number of Lowest-Scoring Species Removed ($k$)')
        ax.set_ylabel('Score $s$')
        ax.set_xlim(0, total_taxa)
        ax.legend(loc='center left', bbox_to_anchor=(0.01, 0.62))
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "combined_score_knee.png"))
        plt.close(fig)

        # Plot 2: Derivative of the sorted score and the cumulative rate deficit it is judged by
        rate_deficit = np.concatenate(([0], np.cumsum(average_rate - score_slope)))
        slope_trend = pd.Series(score_slope).rolling(51, center=True, min_periods=1).median()
        fig, (ax_slope, ax_deficit) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                                                   gridspec_kw={'height_ratios': [3, 2]})
        ax_slope.plot(x[1:], score_slope, color=COLORS['primary'], linewidth=0.6, alpha=0.3,
                      label=r'Slope $\Delta s_k$')
        ax_slope.plot(x[1:], slope_trend, color=COLORS['primary'], label='Rolling median of slope (display only)')
        ax_slope.axhline(average_rate, color=COLORS['secondary'], linewidth=1.5, linestyle=':',
                         label=rf'Average rate $\bar{{m}} = {_sci(average_rate)}$')
        ax_slope.set_yscale('log', nonpositive='mask')
        ax_slope.set_ylabel('Rate of change (log scale)')
        ax_slope.legend(loc='lower center')
        _equation_box(ax_slope, r'$\Delta s_k = s_{k+1} - s_k$' '\n'
                      r'$\bar{m} = \dfrac{s_N - s_1}{N - 1}$', fontsize=12)
        _titles(ax_slope, 'Derivative of Combined Score',
                'Threshold: where the slope overtakes the average rate on balance')

        ax_deficit.fill_between(x, rate_deficit, color=COLORS['primary'], alpha=0.1, linewidth=0)
        ax_deficit.plot(x, rate_deficit, color=COLORS['primary'])
        ax_deficit.plot([k_star], [rate_deficit[knee_idx]], 'o', markersize=9, color=COLORS['primary'],
                        markeredgecolor=COLORS['surface'], markeredgewidth=2, zorder=5)
        ax_deficit.annotate(f'$k^* = {_num(k_star)}$', xy=(k_star, rate_deficit[knee_idx]),
                            xytext=(-14, -16), textcoords='offset points', ha='right', va='top', fontsize=12)
        ax_deficit.set_ylim(bottom=0)
        ax_deficit.set_ylabel('Cumulative rate deficit $D_k$')
        ax_deficit.set_xlabel('Number of Lowest-Scoring Species Removed ($k$)')
        _equation_box(ax_deficit, r'$D_k = \sum_{j<k} \, (\bar{m} - \Delta s_j)$' '\n'
                      r'$k^* = \mathrm{argmax}_k \; D_k$', fontsize=12)
        for ax in (ax_slope, ax_deficit):
            ax.axvline(k_star, color=COLORS['text_secondary'], linewidth=1, linestyle='--')
        ax_deficit.set_xlim(0, total_taxa)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "combined_score_derivative.png"))
        plt.close(fig)

        # Plot 3: Prevalence vs abundance with the decision boundary s = s*
        fig, ax = plt.subplots(figsize=(10, 6))
        removed, retained = taxa_metrics[~keep], taxa_metrics[keep]
        ax.scatter(removed['mean reads'], removed['percent present'], s=10, color=COLORS['removed'],
                   alpha=0.5, linewidths=0, label=f'Removed ({n_removed:,})')
        ax.scatter(retained['mean reads'], retained['percent present'], s=14, color=COLORS['primary'],
                   alpha=0.85, edgecolors=COLORS['surface'], linewidths=0.4, label=f'Retained ({n_kept:,})')
        reads_grid = np.logspace(np.log10(10 ** s_star - 1), np.log10(taxa_metrics['mean reads'].max()), 400)
        ax.plot(reads_grid, 100 * s_star / np.log10(reads_grid + 1), color=COLORS['secondary'],
                label=f'Decision boundary $s = s^* = {s_star:.3f}$')
        ax.set_xscale('log')
        ax.set_ylim(-2, 100)
        _equation_box(ax, r'Keep taxon $i$ if  $p_i \, \log_{10}(\bar{r}_i + 1) > s^*$', x=0.98, y=0.96,
                      ha='right', fontsize=12)
        _titles(ax, 'Species Retained by Combined Score', summary)
        ax.set_xlabel(r'Mean reads $\bar{r}$ (log scale)')
        ax.set_ylabel('Percent of samples present (%)')
        ax.legend(loc='upper left', markerscale=1.5)
        fig.tight_layout()
        fig.savefig(os.path.join(output_dir, "combined_score_retained.png"))
        plt.close(fig)


def analyze_normalized_reads(input_file, save_plots=False, save_filtered=False, output_dir="."):
    df = load_species_table(input_file)
    taxa_metrics = compute_taxa_metrics(df)

    if save_filtered or save_plots:
        os.makedirs(output_dir, exist_ok=True)

    keep, _stats = combined_score_threshold(taxa_metrics, save_plots=save_plots, output_dir=output_dir)

    # Retained taxa with their metrics
    filtered_taxa = taxa_metrics[keep].reset_index().rename(columns={'species': 'Taxon'})

    if save_filtered:
        filtered_taxa.to_csv(os.path.join(output_dir, "filtered_taxa.csv"), index=True)

    return _stats, filtered_taxa


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determine filtering threshold based on species prevalence and read abundance.")
    parser.add_argument("input_file", help="Path to normalized read counts CSV file (species x samples)")
    parser.add_argument("--save-plots", action="store_true", help="Save output plots (default: no plots)")
    parser.add_argument("--save-filtered", action="store_true", help="Save filtered taxa to CSV file (default: no)")
    parser.add_argument("--output-dir", default=".",
                        help="Directory for saved plots and CSV (default: current directory)")
    args = parser.parse_args()

    stats, filtered_taxa = analyze_normalized_reads(args.input_file,
                                                    save_plots=args.save_plots,
                                                    save_filtered=args.save_filtered,
                                                    output_dir=args.output_dir)
    print(f"Combined score threshold: {stats['score_threshold']:.4f}")
    print(f"Species removed: {stats['species_removed']} out of {stats['total_taxa']} "
          f"({stats['species_pct_removed']:.2f}%)")
    print(f"Reads removed: {stats['reads_pct_removed']:.2f}%")
    print(f"Minimum mean reads among retained species: {stats['min_reads_retained']:.1f}")
    print(f"Minimum prevalence among retained species: {stats['min_prevalence_retained']:.2f}%")
