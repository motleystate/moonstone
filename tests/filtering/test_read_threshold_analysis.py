import os
import tempfile
from unittest import TestCase

import numpy as np
import pandas as pd

from moonstone.filtering.read_threshold_analysis import (
    combined_score_threshold, compute_taxa_metrics, load_species_table, rate_threshold_index
)


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
            "Serratia,Serratia_marcescens,615,5,\n"
            "Lactobacillus,Lactobacillus_crispatus,47770,200,300\n")
        table = load_species_table(path)
        self.assertListEqual(table.columns.tolist(), ['s1', 's2'])
        self.assertListEqual(table.loc[('Serratia_marcescens', 615)].tolist(), [15, 0])

    def test_load_species_table_missing_required_column(self):
        path = self._write("genus,species,s1\nSerratia,Serratia_marcescens,10\n")
        with self.assertRaisesRegex(ValueError, 'NCBI_taxonomy_ID'):
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
