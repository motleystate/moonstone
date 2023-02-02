import pandas as pd

from unittest import TestCase

from moonstone.utils.phylogenetic_tree_editing import (
    generate_translation_dictionary,
    replacing_labels,
    adapt_phylogenetic_tree_to_counts_df
)


class TestPhylogeneticTreeAdaptation(TestCase):
    def setUp(self):
        self.count_df = pd.DataFrame(
            [
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillales (order)', 'Lactobacillales (order)', 'Lactobacillales (order)', 186826, 4.3],  # noqa
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_jensenii', 109790, 1.0],  # noqa
                ['Bacteria', 'Firmicutes', 'Bacilli', 'Lactobacillales', 'Lactobacillaceae', 'Lactobacillus', 'Lactobacillus_iners', 147802, 3.5]  # noqa
            ],
            columns=[
                'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species',
                'NCBI_taxonomy_ID', 'SAMPLE_1'
            ]
        )
        self.count_df = self.count_df.set_index(['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'])

    def test_generate_translation_dictionary(self):
        expected_dict = {
            '109790': 'Lactobacillus_jensenii',
            '147802': 'Lactobacillus_iners',
        }
        tested_dict = generate_translation_dictionary(self.count_df['NCBI_taxonomy_ID'])
        self.assertDictEqual(tested_dict, expected_dict)

    def test_replacing_labels(self):
        tree_string = "((('Lactobacillus jensenii, 109790':0.35,\
'Lactobacillus iners, 147802':0.15):0.75,\
'Lactobacillus ruminis CAG:367, 1263085*':1):0.5,\
('Prevotella sp. oral taxon 473, 712469':0.5,\
'Enterococcus lactis, 357441':0.05):1)root;\n"
        tested_string = replacing_labels(
            tree_string,
            {
                '109790': 'Lactobacillus_jensenii',
                '147802': 'Lactobacillus_iners',
                '712469': 'Alloprevotella_Prevotella sp. oral taxon 473',
                '1263085': 'Lactobacillus_ruminis CAG:367'
            }
        )
        expected_string = "((('Lactobacillus_jensenii':0.35,\
'Lactobacillus_iners':0.15):0.75,\
'Lactobacillus_ruminis CAG:367':1):0.5,\
('Alloprevotella_Prevotella sp. oral taxon 473':0.5,\
'Enterococcus lactis, 357441':0.05):1)root;\n"
        self.assertEqual(tested_string, expected_string)

    def test_adapt_phylogenetic_tree_to_counts_df(self):
        tree_string = "((('Lactobacillus jensenii, 109790':0.35,\
'Lactobacillus iners, 147802':0.15):0.75,\
'Lactobacillus ruminis CAG:367, 1263085*':1)root;\n"
        tested_string = adapt_phylogenetic_tree_to_counts_df(
            self.count_df['NCBI_taxonomy_ID'], tree_string
        )
        expected_string = "((('Lactobacillus_jensenii':0.35,\
'Lactobacillus_iners':0.15):0.75,\
'Lactobacillus ruminis CAG:367, 1263085*':1)root;\n"
        self.assertEqual(tested_string, expected_string)
