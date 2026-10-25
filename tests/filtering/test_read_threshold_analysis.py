import contextlib
import io
import os
import tempfile
from unittest import TestCase

import matplotlib
import numpy as np
import pandas as pd

from moonstone.filtering.read_threshold_analysis import (
    analyze_normalized_reads, combined_score_threshold, compute_taxa_metrics, load_species_table, main,
    rate_threshold_index
)

matplotlib.use('Agg')


class TestReadThresholdAnalysis(TestCase):

    def test_compute_taxa_metrics(self):
        tested_object = pd.DataFrame.from_dict(
            {
                'specie_1': [10, 0, 0, 2],
                'specie_2': [4, 4, 4, 4],
            },
            orient='index', columns=['1', '2', '3', '4'])
        metrics = compute_taxa_metrics(tested_object)
        self.assertListEqual(metrics['percent present'].tolist(), [50.0, 100.0])
        self.assertListEqual(metrics['mean reads'].tolist(), [3.0, 4.0])
        self.assertListEqual(metrics['total reads'].tolist(), [12, 16])

    def test_rate_threshold_index(self):
        # Flat then sharply rising: threshold is the last point before the rise
        y = [0, 0, 0, 0, 0, 0, 0, 0, 1, 10]
        idx, slope, average_rate = rate_threshold_index(y)
        self.assertEqual(idx, 8)
        self.assertEqual(average_rate, 10 / 9)
        self.assertListEqual(slope.tolist(), [0, 0, 0, 0, 0, 0, 0, 1, 9])

    def test_rate_threshold_index_ignores_early_noisy_step(self):
        # One early step above the average rate should not set the threshold
        y = [0, 0, 3, 3, 3, 3, 3, 3, 4, 20]
        idx, slope, average_rate = rate_threshold_index(y)
        self.assertGreater(slope[1], average_rate)
        self.assertEqual(idx, 8)

    def test_combined_score_requires_prevalence_and_abundance(self):
        taxa_metrics = pd.DataFrame(
            {
                'percent present': [90, 85, 80, 2, 90, 1, 1, 1, 1, 1],
                'mean reads': [500, 400, 300, 1000, 0.1, 0.1, 0.2, 0.3, 0.1, 0.2],
            },
            index=[f'specie_{i}' for i in range(10)])
        keep, stats = combined_score_threshold(taxa_metrics)
        self.assertEqual(int(keep.sum()), 3)
        # Prevalent and abundant taxa are kept
        self.assertTrue(keep[['specie_0', 'specie_1', 'specie_2']].all())
        # A single-sample spike and a prevalent trace taxon are not
        self.assertFalse(keep['specie_3'])
        self.assertFalse(keep['specie_4'])
        self.assertEqual(stats['species_removed'], int((~keep).sum()))
        self.assertTrue(np.isclose(stats['min_reads_retained'], 300))


class TestLoadSpeciesTable(TestCase):

    def _write(self, text):
        handle = tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False)
        handle.write(text)
        handle.close()
        self.addCleanup(os.remove, handle.name)
        return handle.name

    def test_load_species_table(self):
        path = self._write(
            "genus,species,NCBI_taxonomy_ID,s1,s2\n"
            "Serratia,Serratia_marcescens,615,10,0\n"
            "Serratia,Serratia_marcescens,615,5,0\n"
            "Lactobacillus,Lactobacillus_crispatus,47770,200,300\n")
        table = load_species_table(path)
        self.assertListEqual(table.columns.tolist(), ['s1', 's2'])
        self.assertListEqual(table.loc[('Serratia_marcescens', 615)].tolist(), [15, 0])

    def test_load_species_table_blank_count_raises(self):
        path = self._write("species,NCBI_taxonomy_ID,s1,s2\nSerratia_marcescens,615,5,\n")
        with self.assertRaisesRegex(ValueError, '1 blank count cell.*s2'):
            load_species_table(path)

    def test_load_species_table_keeps_rows_without_id(self):
        path = self._write("species,NCBI_taxonomy_ID,s1\nSerratia_marcescens,615,5\nUnknown_sp,,500\n")
        with self.assertWarnsRegex(UserWarning, '1 row\\(s\\) with a blank species or NCBI_taxonomy_ID'):
            table = load_species_table(path)
        self.assertEqual(table.shape[0], 2)
        self.assertEqual(table['s1'].sum(), 505)

    def test_load_species_table_missing_required_column(self):
        path = self._write("genus,species,s1\nSerratia,Serratia_marcescens,10\n")
        with self.assertRaisesRegex(ValueError, 'NCBI_taxonomy_ID'):
            load_species_table(path)

    def test_load_species_table_no_sample_columns(self):
        path = self._write("genus,species,NCBI_taxonomy_ID\nSerratia,Serratia_marcescens,615\n")
        with self.assertRaisesRegex(ValueError, 'no sample columns'):
            load_species_table(path)

    def test_load_species_table_non_numeric_column(self):
        path = self._write("species,NCBI_taxonomy_ID,subspecies,s1\nSerratia_marcescens,615,abc,10\n")
        with self.assertRaisesRegex(ValueError, 'subspecies'):
            load_species_table(path)

    def test_load_species_table_negative_values(self):
        path = self._write("species,NCBI_taxonomy_ID,s1\nSerratia_marcescens,615,-1\n")
        with self.assertRaisesRegex(ValueError, 'negative'):
            load_species_table(path)

    def test_load_species_table_warns_on_relative_abundances(self):
        path = self._write("species,NCBI_taxonomy_ID,s1,s2\nA,1,0.4,0.9\nB,2,0.6,0.1\n")
        with self.assertWarnsRegex(UserWarning, 'relative abundances'):
            load_species_table(path)


class TestAnalyzeNormalizedReads(TestCase):
    """End-to-end: CSV in, retained taxa, plots and filtered CSV out."""

    def setUp(self):
        self.tmp_dir = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp_dir.cleanup)
        samples = [f's{i}' for i in range(1, 11)]
        rows = {
            # Prevalent and abundant: kept
            'Abundant_A': [5000] * 10,
            'Abundant_B': [2000] * 9 + [0],
            'Abundant_C': [800] * 8 + [0] * 2,
            # Abundant in a single sample only: removed
            'Spike': [50000] + [0] * 9,
        }
        # Single reads in single samples: removed
        rows.update({f'Trace_{i}': [0] * i + [1] + [0] * (9 - i) for i in range(6)})
        df = pd.DataFrame.from_dict(rows, orient='index', columns=samples)
        df.insert(0, 'NCBI_taxonomy_ID', range(1, len(df) + 1))
        df.insert(0, 'species', df.index)
        df.insert(0, 'genus', 'Genus')
        self.input_file = os.path.join(self.tmp_dir.name, 'counts.csv')
        df.to_csv(self.input_file, index=False)
        self.output_dir = os.path.join(self.tmp_dir.name, 'results')

    def test_analyze_normalized_reads(self):
        stats, filtered_taxa = analyze_normalized_reads(self.input_file, save_plots=True, save_filtered=True,
                                                        output_dir=self.output_dir)
        self.assertListEqual(sorted(filtered_taxa['Taxon']), ['Abundant_A', 'Abundant_B', 'Abundant_C'])
        self.assertEqual(stats['total_taxa'], 10)
        self.assertEqual(stats['species_removed'], 7)
        self.assertEqual(stats['min_reads_retained'], 640)
        self.assertEqual(stats['min_prevalence_retained'], 80)
        self.assertListEqual(sorted(os.listdir(self.output_dir)),
                             ['combined_score_derivative.png', 'combined_score_knee.png',
                              'combined_score_retained.png', 'filtered_taxa.csv'])
        saved = pd.read_csv(os.path.join(self.output_dir, 'filtered_taxa.csv'))
        self.assertListEqual(sorted(saved['Taxon']), ['Abundant_A', 'Abundant_B', 'Abundant_C'])

    def test_analyze_normalized_reads_writes_nothing_by_default(self):
        analyze_normalized_reads(self.input_file, output_dir=self.output_dir)
        self.assertFalse(os.path.exists(self.output_dir))

    def test_command_line(self):
        stdout = io.StringIO()
        with contextlib.redirect_stdout(stdout):
            main([self.input_file, '--save-filtered', '--output-dir', self.output_dir])
        self.assertIn('Species removed: 7 out of 10', stdout.getvalue())
        self.assertTrue(os.path.exists(os.path.join(self.output_dir, 'filtered_taxa.csv')))
