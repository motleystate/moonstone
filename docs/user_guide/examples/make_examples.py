"""
Write the synthetic example files used in the user guide (docs/user_guide/correlation.rst).

All values are simulated: 30 patients sampled at inclusion and at delivery. Gardnerella vaginalis, Fannyhessea
vaginae and Prevotella bivia are higher, and Lactobacillus crispatus lower, in patients with bacterial vaginosis (bv);
the resistance gene tetM follows Gardnerella vaginalis. Run from this folder: python make_examples.py
"""
import numpy as np
import pandas as pd

TAXONOMY = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
# taxonomy, mean log abundance, shift with bv
TAXA = [
    (('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus',
      'Lactobacillus crispatus'), 8.0, -2.5),
    (('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus',
      'Lactobacillus iners'), 6.5, 0.0),
    (('Bacteria', 'Actinobacteria', 'Actinomycetia', 'Bifidobacteriales', 'Bifidobacteriaceae', 'Gardnerella',
      'Gardnerella vaginalis'), 4.5, 2.5),
    (('Bacteria', 'Actinobacteria', 'Coriobacteriia', 'Coriobacteriales', 'Atopobiaceae', 'Fannyhessea',
      'Fannyhessea vaginae'), 3.0, 2.0),
    (('Bacteria', 'Bacteroidetes', 'Bacteroidia', 'Bacteroidales', 'Prevotellaceae', 'Prevotella',
      'Prevotella bivia'), 3.5, 1.5),
    (('Bacteria', 'Fusobacteria', 'Fusobacteriia', 'Fusobacteriales', 'Leptotrichiaceae', 'Sneathia',
      'Sneathia amnii'), 1.0, 1.0),
    (('Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Streptococcaceae', 'Streptococcus',
      'Streptococcus agalactiae'), 1.5, 0.0),
    (('Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacterales', 'Enterobacteriaceae', 'Escherichia',
      'Escherichia coli'), 1.0, 0.0),
]


def main(n_patients=30, seed=2026):
    rng = np.random.default_rng(seed)
    n = 2 * n_patients
    samples = [f'S{i:03d}' for i in range(1, n + 1)]
    patient = np.repeat(np.arange(n_patients), 2)
    bv = rng.permutation([True] * 12 + [False] * (n_patients - 12))[patient]
    metadata = pd.DataFrame({
        'sample': samples,
        'patient': [f'P{p + 1:02d}' for p in patient],
        'Type': ['inclusion', 'delivery'] * n_patients,
        'disease': np.where(bv, 'bv', 'healthy'),
        'age': rng.integers(20, 41, n_patients)[patient],
        'antibiotics': rng.choice(['yes', 'no'], n_patients, p=[0.3, 0.7])[patient],
    })
    log_abundance = np.array([mean + shift * bv + rng.normal(0, 1, n) for _, mean, shift in TAXA])
    proportions = np.exp(log_abundance) / np.exp(log_abundance).sum(axis=0)
    counts = np.stack([rng.multinomial(rng.integers(15000, 40000), proportions[:, j]) for j in range(n)]).T
    species = pd.concat([pd.DataFrame([taxonomy for taxonomy, _, _ in TAXA], columns=TAXONOMY),
                         pd.DataFrame(counts, columns=samples)], axis=1)
    species.to_csv('species_counts.tsv', sep='\t', index=False)

    genes = pd.DataFrame(index=pd.Index(['tetM', 'ermB', 'lnuA', 'blaTEM'], name='gene'), columns=samples, dtype=float)
    genes.loc['tetM'] = np.round(counts[2] / counts.sum(axis=0) * 200 * rng.lognormal(0, 0.3, n), 1)
    for gene, mean in [('ermB', 3.0), ('lnuA', 1.5), ('blaTEM', 0.8)]:
        values = np.round(rng.lognormal(mean, 1, n), 1)
        values[rng.random(n) < 0.25] = 0
        genes.loc[gene] = values
    genes.to_csv('resistance_genes.csv')
    metadata.to_csv('metadata.csv', index=False)


if __name__ == '__main__':
    main()
