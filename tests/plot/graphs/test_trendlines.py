from unittest import TestCase

import pandas as pd

from moonstone.plot.graphs.trendlines import ScatterTrendlines


class TestScatterTrendlines(TestCase):

    def setUp(self):
        sp1 = [
            1.93, 0.46, 1.82, 1.27, 0.39, 1.43, 0.86, 0.64, 1.86, 2.23, 0.35, 1.11, 2.84, 2.59, 0.86, 0.89, 2.29,
            2.73, 0.1, 2.47, 2.38, 0.85, 0.11, 0.57, 0.04, 2.79, 1.32, 1.74, 0.97, 2.55, 2.61, 2.29, 2.73, 0.83, 1.99,
            2.15, 0.45, 1.86, 0.19, 2.23, 1.94, 0.16, 2.12, 2.1, 0.49, 0.09, 2.46, 1.29, 2.27, 1.71, 1.08, 2.68, 2.46,
            1.07, 0.04, 2.93, 2.82, 0.82, 0.01, 0.6, 1.64, 2.42, 0.87, 0.94, 2.92, 0.32, 2.83, 1.85, 1.53, 0.61]
        sp2 = [
            74.48, 63.57, 76.55, 30.21, 23.32, 66.4, 58.17, 56.32, 35.61, 71.73, 18.94, 66.58, 85.25, 86.56, 26.56,
            23.62, 36.68, 77.23, 47.07, 92.92, 75.14, 65.14, 9.84, 62.44, 44.44, 38.71, 27.16, 36.79, 63.28, 47.42,
            36.51, 73.28, 81.97, 62.44, 64.38, 68.21, 58.81, 33.09, 24.01, 34.02, 73.77, 18.23, 39.68, 33.33, 62.37,
            19.21, 77.23, 74.69, 80.31, 24.07, 24.84, 46.0, 69.46, 62.2, 24.27, 86.64, 85.0, 25.09, 56.48, 25.72,
            80.58, 75.11, 27.22, 73.6, 48.52, 22.24, 42.2, 83.64, 73.16, 55.14]
        self.samples = [
            'A11', 'A12', 'A16', 'B3', 'B11', 'A27', 'A23', 'A20', 'B4', 'A8', 'B17', 'A30', 'A3', 'A29', 'B26', 'B20',
            'B21', 'A14', 'A31', 'A39', 'A37', 'A0', 'B15', 'A6', 'A4', 'B7', 'B27', 'B24', 'A38', 'B9', 'B25', 'A2',
            'A18', 'A28', 'A1', 'A21', 'A7', 'B5', 'B19', 'B28', 'A15', 'B22', 'B14', 'B16', 'A19', 'B23', 'A34',
            'A33', 'A24', 'B6', 'B13', 'B1', 'A17', 'A25', 'B0', 'A36', 'A10', 'B29', 'A9', 'B2', 'A5', 'A35', 'B12',
            'A13', 'B18', 'B8', 'B10', 'A22', 'A32', 'A26']
        self.test_df = pd.DataFrame(
            [sp1, sp2], columns=self.samples, index=["species1", "species2"]
        ).T
        self.expected_group_series = self.test_df.index.to_series().map(lambda x: 0 if x[0] == "A" else 1)
        self.expected_group_series.name = "group"
        self.ins = ScatterTrendlines(self.test_df)

    def test_define_n_trendlines(self):
        output = self.ins.define_n_trendlines(
            "species1", "species2", 2, False, False, nb_iter=100, outliers='keep', show=False, v=False,
            random_state=2
        )
        # output =
        #   * group_series
        pd.testing.assert_series_equal(output["group_series"], self.expected_group_series)
        #   * distance_dataframe
        dist = [
            0.0429662986986339, 0.4858022223403532, 0.34614862434873656, 0.129570767158105, 0.22689243838279075,
            0.2151744165543507, 0.41765584418798435, 0.371702869050711, 0.152354574660951, 0.5130020126313728,
            0.22734525572914377, 0.12026243077273509, 0.14457602438112896, 0.5160403023926496, 0.1253028230974578,
            0.23608305756485357, 0.15422586337017322, 0.496240616994477, 0.6994817453033855, 1.2305420537388476,
            0.34331173249677155, 0.24439868870621806, 1.0151871443508251, 0.2705642999612276, 0.8858021281455863,
            0.42208955503572365, 0.2640934628926417, 0.4046709936659474, 0.0490926214491288, 0.7987051079403361,
            0.4913570438692193, 0.4277241356749775, 0.05277628560678041, 0.011704700846672293, 0.9617050776501762,
            0.7626766280985767, 0.05042287186373525, 0.1318545275303066, 0.5034355624661815, 0.39460717058259076,
            0.033415931561183956, 0.11863227556506659, 0.3530336762826934, 0.3432557521683666, 0.3436643560736769,
            0.061446874200724706, 0.22742487945243897, 0.6998063119797779, 0.24990000819591857, 1.0000996892888445,
            0.28727752895212877, 0.5093849571358114, 0.9543695745753105, 0.2496965116097618, 0.681801697665927,
            0.1850162259238461, 0.1410988322810478, 0.0007410268174285378, 0.2705031120119411, 0.28890761039546287,
            0.9023973555987914, 0.38594302248130835, 0.18980234208052188, 0.9462930722976726, 0.5551253011469858,
            0.1746419252462042, 0.06822729780879641, 0.9796055486314305, 0.3177154397152984, 0.45223274723608536]
        expected_distance_dataframe = pd.DataFrame([dist], columns=self.samples, index=["distance"]).T
        pd.testing.assert_frame_equal(output["distance_dataframe"], expected_distance_dataframe)
        #   * trendlines
        self.assertEqual(len(output["trendlines"]), 2)
        # checking trendline of the first group
        self.assertEqual(output["trendlines"][0][1], 0)
        x = [
            1.64, 0.46, 0.97, 2.73, 2.82, 2.15, 0.57, 2.42, 1.93, 2.27, 0.49, 1.85, 0.45, 2.23, 2.59, 0.1, 2.73, 1.53,
            0.61, 2.84, 1.82, 2.46, 2.47, 0.01, 0.86, 0.04, 1.43, 2.93, 0.83, 0.94, 1.07, 2.46, 2.38, 1.29, 0.85,
            1.11, 2.29, 1.99, 0.64, 1.94]
        y = [
            70.93466177820105, 58.37746844579898, 63.80472997082021, 82.53410307677584, 83.4918551106031,
            76.36192330322227, 59.548054264921205, 79.2351794047041, 74.02075166497782, 77.63892601499198,
            58.69671912374141, 73.16941652379802, 58.27105155315151, 77.21325844440207, 81.04426657971118,
            54.54646031048988, 82.53410307677584, 69.76407595907882, 59.97372183551111, 83.70468889589806,
            72.8501658458556, 79.660846975294, 79.76726386794148, 53.5887082766626, 62.63414415169799,
            53.90795895460503, 68.69990703260407, 84.66244092972534, 62.31489347375556, 63.485479292877784,
            64.86889889729497, 79.660846975294, 78.8095118341142, 67.21007053553942, 62.52772725905051,
            65.29456646788486, 77.85175980028693, 74.65925302086268, 60.292972513453535, 74.1271685576253]
        samples_group0 = [
            'A5', 'A12', 'A38', 'A18', 'A10', 'A21', 'A6', 'A35', 'A11', 'A24', 'A19', 'A22', 'A7', 'A8', 'A29', 'A31',
            'A14', 'A32', 'A26', 'A3', 'A16', 'A34', 'A39', 'A9', 'A23', 'A4', 'A27', 'A36', 'A28', 'A13', 'A25',
            'A17', 'A37', 'A33', 'A0', 'A30', 'A2', 'A1', 'A20', 'A15']
        expected_group0_trendlines_coor = pd.DataFrame(
            [x, y], columns=samples_group0, index=["x", "y"]
        ).T
        pd.testing.assert_frame_equal(output["trendlines"][0][0], expected_group0_trendlines_coor)
        #   * figure
        # checking that it exists
        self.assertEqual(output["figure"].data[0]["type"], "scatter")

    def test_define_n_trendlines_trim_outliers(self):
        ins = ScatterTrendlines(
            pd.concat([
                self.test_df,
                pd.DataFrame([[0, 100]], columns=["species1", "species2"], index=["outliers1"])
                ]))
        output = ins.define_n_trendlines(
            "species1", "species2", 2, False, False, nb_iter=100, outliers='trim', show=False, v=False,
            random_state=2
        )
        expected_group_series = pd.concat([self.expected_group_series, pd.Series([0], index=["outliers1"])])
        expected_group_series.name = "group"
        pd.testing.assert_series_equal(output["group_series"], expected_group_series)

    def test_define_n_trendlines_exception1(self):
        # Exception 1 = case where less than 2 points in one of the chunk_lists
        # which makes it impossible to compute a (trend)line
        output = self.ins.define_n_trendlines(
            "species1", "species2", 3, False, False, nb_iter=100, outliers='keep', show=False, v=False,
            random_state=19
        )
        pd.testing.assert_series_equal(output["group_series"], self.expected_group_series)
        self.assertEqual(len(output["trendlines"]), 2)

    def test_define_n_trendlines_exception2(self):
        # Exception 2 = case where n=3 trendlines where created
        # but during the assignment of the datapoints to one of the trendlines,
        # one of the trendline didn't get assigned any datapoints
        output = self.ins.define_n_trendlines(
            "species1", "species2", 3, False, False, nb_iter=100, outliers='keep', show=False, v=False,
            random_state=5
        )
        pd.testing.assert_series_equal(output["group_series"], self.expected_group_series)
        self.assertEqual(len(output["trendlines"]), 2)

    def test_bootstraps_define_n_trendlines(self):
        output = self.ins.bootstraps_define_n_trendlines(
            "species1", "species2", 2, False, False, nb_iter=100, nb_bootstraps=5, outliers="keep", show=False,
            random_state=11
        )
        pd.testing.assert_series_equal(output["group_series"], self.expected_group_series)

    def test_plot_one_graph(self):
        # no bootstraps -> calls define_n_trendlines()
        fig = self.ins.plot_one_graph(
            "species1", "species2", 2, log=False, nb_iter=100, nb_bootstraps=1, outliers="keep", show=False,
            random_state=2
        )
        # just checking that it runs with errors and returns something
        self.assertEqual(fig.data[0]["type"], "scatter")

        # 2 bootstraps -> calls bootstraps_define_n_trendlines()
        fig = self.ins.plot_one_graph(
            "species1", "species2", 2, log=False, nb_iter=100, nb_bootstraps=2, outliers="keep", show=False,
            random_state=2
        )
        self.assertEqual(fig.data[0]["type"], "scatter")

    def test_color_palette(self):
        # simplest case
        tested_object = self.ins._color_palette(5)
        self.assertEqual(len(tested_object), 5)

        tested_object = self.ins._color_palette(10)
        self.assertEqual(len(tested_object), 10)

        # n > 10 -> second palette
        tested_object = self.ins._color_palette(15)
        self.assertEqual(len(tested_object), 15)

        # n > 26 -> second palette repeated over and over again
        tested_object = self.ins._color_palette(123)
        self.assertEqual(len(tested_object), 123)

#    def test_define_n_trendlines_log(self):
#       @TODO: write test with x and y as log
