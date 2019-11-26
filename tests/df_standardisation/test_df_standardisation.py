from unittest import TestCase

import pandas as pd

from moonstone.dataframe_standardisation.df_standardisation import (
    Df_Standardisation
)


class TestDf_Standardisation(TestCase):
    """
    expected data filepath = exdfp
    tested data filepath = tdfp
    """

    def test_import_df(self):
        tdfp = '/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/raw_data.csv'
        tested_object = Df_Standardisation(tdfp)
        exdfp = '/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/test_1_import_data.csv'
        expected_object = pd.read_csv(exdfp, dtype=str)
        expected_object.index += 1
        pd.testing.assert_frame_equal(tested_object.import_df(tdfp), expected_object)

    def test_delete_lowercase_uppercase(self):
        test_string = "Salamandra"
        tested_object = Df_Standardisation(test_string)
        expected_object = 'Salamandra'
        self.assertEqual(tested_object.delete_lowercase(test_string), expected_object)

    def test_delete_lowercase_lowercase(self):
        test_string = "salamandra"
        tested_object = Df_Standardisation(test_string)
        expected_object = None
        self.assertEqual(tested_object.delete_lowercase(test_string), expected_object)

    def test_spliting_into_taxa_columns(self):
        filepath = '/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/test_1_import_data.csv'
        tdfp = pd.read_csv(filepath)
        tested_object = Df_Standardisation(tdfp)
        exdfp = '/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/test_2_spliting_taxa_columns.csv'
        exdfp = [
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Bacteroidetes',	None, 'D', 2, "", 'Bacteroidia',
             "D", 3, "", 'Bacteroidales', 'D', 4, "", 'Tannerellaceae', 'D', 5, "", 'Macellibacteroides'],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Proteobacteria',	None, 'D', 2, "", 'Gammaproteobacteria',
             "D", 3, "", 'Betaproteobacteriales', 'D', 4, "", 'Chitinibacteraceae', 'D', 5, "", 'Deefgea'],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Actinobacteria',	None, 'D', 2, "", 'Acidimicrobiia',
             "D", 3, "", 'Microtrichales', 'D', 4, "", 'Microtrichaceae', 'D', 5, "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Deinococcus',	'Thermus', 'D', 2, "", 'Deinococci',
             "D", 3, "", 'Deinococcales', 'D', 4, "", 'Deinococcaceae', 'D', 5, "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Proteobacteria',	None, 'D', 2, "", 'Alphaproteobacteria',
             "D", 3, "", 'Caulobacterales', 'D', 4, "", 'Hyphomonadaceae', "", "", "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Patescibacteria',	None, 'D', 2, "", 'Saccharimonadia',
             "D", 3, "", 'Saccharimonadales', 'D', 4, "", None, 'D', 5, "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Proteobacteria',	None, 'D', 2, "", 'Alphaproteobacteria',
             "D", 3, "", 'Rhodobacterales', 'D', 4, "", 'Rhodobacteraceae', 'D', 5, "", 'Cereibacter'],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Fibrobacteres',	None, 'D', 2, "", 'Fibrobacteria',
             "D", 3, "", 'Fibrobacterales', 'D', 4, "", 'Fibrobacteraceae', 'D', 5, "", None],
        ]
        expected_object = pd.DataFrame(exdfp, columns=[0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 3,
                                                       0, 1, 2, 3, 0, 1, 2, 3], dtype=str)
        pd.testing.assert_frame_equal(tested_object.spliting_into_taxa_columns(tdfp), expected_object)

    def test_naming_taxa_columns(self):
        tdfp = [
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Bacteroidetes',	None, 'D', 2, "", 'Bacteroidia',
             "D", 3, "", 'Bacteroidales', 'D', 4, "", 'Tannerellaceae', 'D', 5, "", 'Macellibacteroides'],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Proteobacteria',	None, 'D', 2, "", 'Gammaproteobacteria',
             "D", 3, "", 'Betaproteobacteriales', 'D', 4, "", 'Chitinibacteraceae', 'D', 5, "", 'Deefgea'],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Actinobacteria',	None, 'D', 2, "", 'Acidimicrobiia',
             "D", 3, "", 'Microtrichales', 'D', 4, "", 'Microtrichaceae', 'D', 5, "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Deinococcus',	'Thermus', 'D', 2, "", 'Deinococci',
             "D", 3, "", 'Deinococcales', 'D', 4, "", 'Deinococcaceae', 'D', 5, "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Proteobacteria',	None, 'D', 2, "", 'Alphaproteobacteria',
             "D", 3, "", 'Caulobacterales', 'D', 4, "", 'Hyphomonadaceae', "", "", "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Patescibacteria',	None, 'D', 2, "", 'Saccharimonadia',
             "D", 3, "", 'Saccharimonadales', 'D', 4, "", None, 'D', 5, "", None],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Proteobacteria',	None, 'D', 2, "", 'Alphaproteobacteria',
             "D", 3, "", 'Rhodobacterales', 'D', 4, "", 'Rhodobacteraceae', 'D', 5, "", 'Cereibacter'],
            ['D', 0, "", 'Bacteria', 'D', 1, "", 'Fibrobacteres',	None, 'D', 2, "", 'Fibrobacteria',
             "D", 3, "", 'Fibrobacterales', 'D', 4, "", 'Fibrobacteraceae', 'D', 5, "", None],
        ]
        tdfp = pd.DataFrame(tdfp, columns=[0, 1, 2, 3, 0, 1, 2, 3, 4, 0, 1, 2, 3, 0, 1, 2, 3,
                                           0, 1, 2, 3, 0, 1, 2, 3], dtype=str)
        tested_object = Df_Standardisation(tdfp)
        exdfp = "/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/test_3_taxa_data.csv"
        expected_object = pd.read_csv(exdfp, dtype=str)
        pd.testing.assert_frame_equal(tested_object.naming_taxa_columns(tdfp), expected_object)

    def test_filling_missing_taxa_values(self):
        filepath = "/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/test_3_taxa_data.csv"
        tdfp = pd.read_csv(filepath, dtype=str)
        tested_object = Df_Standardisation(filepath)
        exdfp = "/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/test_4_taxa_df_completed.csv"
        expected_object = pd.read_csv(exdfp, dtype=str)
        pd.testing.assert_frame_equal(tested_object.filling_missing_taxa_values(tdfp), expected_object)

    def test_standard_taxa_df(self):
        tdfp = '/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/raw_data.csv'
        tested_object = Df_Standardisation(tdfp)
        print(tested_object.standard_taxa_df)
        exdfp = '/home/mariela/Pasteur/Gitlab/data-analysis/tests/df_standardisation/test_5_final_df_with_taxa.csv'
        expected_object = pd.read_csv(exdfp, index_col=['kingdom', 'phylum', 'class', 'order', 'family', 'genus'],
                                      dtype=str)
        pd.testing.assert_frame_equal(tested_object.standard_taxa_df, expected_object)
