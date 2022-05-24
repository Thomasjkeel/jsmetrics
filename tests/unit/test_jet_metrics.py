# -*- coding: utf-8 -*-

"""
    Tests for jet metrics (jetstream_metrics.py & jetstream_metrics_utils)
    Includes:
        – tests for each metric defined and expected outputs
        – test for each utility function used in the metric
          (e.g. core algorithms)
"""

# imports
import unittest
from parameterized import parameterized
import xarray as xr
import numpy as np
from jsmetrics import (
    general_utils,
    jetstream_metrics,
    jetstream_metrics_utils,
    jetstream_metrics_dict,
)
from . import (
    set_up_test_uv_data,
    set_up_test_u_data,
    set_up_test_zg_data,
    set_up_nan_dataset,
)


# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


MAX_VARIABLES = 7


class TestJetStreamMetricDict(unittest.TestCase):
    def setUp(self):
        self.metric_dict = jetstream_metrics_dict.JETSTREAM_METRIC_DICT

    def test_metric_dict_keys(self):
        for metric_name in self.metric_dict.keys():
            self.assertIsInstance(metric_name, str)

    def test_metric_dict_values(self):
        for metric in self.metric_dict.values():
            self.assertIsInstance(metric, dict)
            self.assertEqual(len(metric.keys()), MAX_VARIABLES)
            self.assertListEqual(
                list(metric.keys()),
                [
                    "variables",
                    "coords",
                    "plev_units",
                    "metric",
                    "name",
                    "description",
                    "doi",
                ],
            )

    def test_variables(self):
        for metric in self.metric_dict.values():
            self.assertIsInstance(metric["variables"], list)
            self.assertGreaterEqual(len(metric["variables"]), 0)
            self.assertLessEqual(len(metric["variables"]), MAX_VARIABLES)

    def test_metric_coords(self):
        for metric in self.metric_dict.values():
            self.assertIsInstance(metric["coords"], dict)
            for coord in metric["coords"].keys():
                self.assertIsInstance(coord, str)

    @parameterized.expand(
        [
            ("plev", 0, 100000),
            ("lat", -91, 91),
            ("lon", -1, 361),
        ]
    )
    def test_each_coord(self, coord, min_value, max_value):
        for metric in self.metric_dict.values():
            if coord in metric["coords"].keys():
                self.assertEqual(len(metric["coords"][coord]), 2)
                self.assertGreaterEqual(
                    min(metric["coords"][coord]), min_value
                )
                self.assertLessEqual(max(metric["coords"][coord]), max_value)

    def test_funcs(self):
        for metric in self.metric_dict.values():
            self.assertTrue(callable(metric["metric"]))


class TestKoch2006(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jetstream_metrics.koch_et_al_2006
        result = tested_func(self.data, ws_threshold=8)
        # check an exact value. Is this necessary?
        self.assertIsInstance(result, xr.Dataset)
        self.assertEqual(
            float(result["weighted_average_ws"].max()), 8.775158882141113
        )
        new_data = self.data.isel(plev=0)
        self.assertRaises(ValueError, lambda: tested_func(new_data))

    def test_get_all_hPa_list(self):
        tested_func = general_utils.get_all_hPa_list
        # make sure it returns an array
        self.assertIsInstance(tested_func(self.data), (np.ndarray))
        # make sure it takes errors wrong types
        # self.assertRaises(TypeError,\
        # lambda: jetstream_metrics_utils.get_all_plev_hPa(['plev']))
        new_data = self.data.rename({"plev": "pl"})
        self.assertRaises(KeyError, lambda: tested_func(new_data))
        new_data2 = self.data.copy()
        new_data2["plev"] = self.data.plev.assign_attrs(units="mbar")
        self.assertRaises(ValueError, lambda: tested_func(new_data2))

    def test_sum_weighted_ws(self):
        tested_func = jetstream_metrics_utils.get_sum_weighted_ws
        self.assertRaises(TypeError, lambda: tested_func(self.data, 1))
        new_data = self.data.rename({"plev": "pl"})
        self.assertRaises(KeyError, lambda: tested_func(new_data, [10, 20]))
        sum_weighted = tested_func(self.data, [0, 100])
        self.assertIsInstance(sum_weighted, xr.DataArray)
        self.assertGreater(sum_weighted.max(), 0)

    def test_weighted_average_ws(self):
        tested_func = jetstream_metrics_utils.get_weighted_average_ws

        self.assertRaises(TypeError, lambda: tested_func(self.data, 1))
        sum_weighted = jetstream_metrics_utils.get_sum_weighted_ws(
            self.data, [0, 100]
        )
        weighted_av = tested_func(sum_weighted, np.array([0, 100]))
        self.assertGreater(weighted_av.min(), 0)
        self.assertGreaterEqual(weighted_av.max(), 0)

        weighted_av_threshold = weighted_av.where(weighted_av >= 30)
        self.assertGreater(weighted_av_threshold.max(), 0)
        self.assertGreaterEqual(weighted_av_threshold.min(), 0)


class TestArcherCaldeira2008(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        result = jetstream_metrics.archer_caldeira_2008(self.data)
        for col in [
            "mass_weighted_average_ws",
            "mass_flux_weighted_pressure",
            "mass_flux_weighted_latitude",
        ]:
            self.assertIn(col, result)
            # self.assertEqual(10, result[col].max())

        self.assertEqual(
            float(result["mass_weighted_average_ws"].max()), 23.904821395874023
        )

    def test_get_mass_weighted_average_windspeed(self):
        tested_func = jetstream_metrics.archer_caldeira_2008
        test_data = self.data.isel(plev=1)
        # should fail because needs two plevs
        self.assertRaises(ValueError, lambda: tested_func(test_data))


class TestSchiemann2009(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jetstream_metrics.schiemann_et_al_2009
        sub_data = self.data.isel(time=slice(0, 1), lon=slice(0, 30))
        result = tested_func(sub_data)
        self.assertTrue(result)


class TestWoolings2010(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()
        # self.data = make_fake_seasonal_data(self.data)
        self.test_sig = self._get_test_sig()

    @staticmethod
    def _get_test_sig():
        # Seed the random number generator
        np.random.seed(42)
        time_step = 0.25
        period = 5.0
        time_vec = np.arange(0, 5, time_step)
        test_sig = np.sin(
            2 * np.pi / period * time_vec
        ) + 0.5 * np.random.randn(time_vec.size)
        return test_sig

    def test_metric(self):
        result = jetstream_metrics.woolings_et_al_2010(
            self.data, filter_freq=1, window_size=2
        )
        self.assertIsInstance(result, xr.Dataset)
        self.assertEqual(result["ff_max_lats"][0], 36.25)
        self.assertEqual(result["ff_max_ws"][0], 43.365413665771484)

    def test_get_zonal_mean(self):
        tested_func = jetstream_metrics_utils.get_zonal_mean
        new_data = self.data.rename({"lon": "ln"})
        self.assertRaises(KeyError, lambda: tested_func(new_data))
        self.assertIsInstance(tested_func(self.data), xr.Dataset)

    def test_apply_lancoz_filter(self):
        tested_func = jetstream_metrics_utils.apply_lanczos_filter
        test_ua = self.data["ua"]
        self.assertRaises(AssertionError, lambda: tested_func(self.data, 2, 4))
        self.assertRaises(AssertionError, lambda: tested_func(test_ua, -2, 1))
        self.assertRaises(ValueError, lambda: tested_func(test_ua, 2, -1))
        self.assertRaises(ValueError, lambda: tested_func(test_ua, 2, 1))
        self.assertRaises(
            ValueError,
            lambda: tested_func(test_ua, test_ua["time"].count() + 2, 1),
        )
        self.assertRaises(
            ValueError,
            lambda: tested_func(test_ua, 2, test_ua["time"].count() + 1),
        )
        self.assertEqual(
            float(tested_func(test_ua, 2, 4).max()), 99.514892578125
        )

    def test_get_latitude_and_speed_where_max_ws(self):
        tested_func = (
            jetstream_metrics_utils.get_latitude_and_speed_where_max_ws
        )
        self.assertRaises(AttributeError, lambda: tested_func(["lol"]))
        tested_data = self.data["ua"].isel(plev=0, lon=0, time=0)
        self.assertEqual(tested_func(tested_data)[0], 81.25)
        self.assertEqual(tested_func(tested_data)[1], -9.87109375)
        self.assertRaises(
            KeyError, lambda: tested_func(tested_data.rename({"lat": "lt"}))
        )
        nan_dataset = set_up_nan_dataset()
        self.assertEqual(tested_func(nan_dataset), (np.nan, np.nan))

    def test_apply_fourier_filter(self):
        tested_func = jetstream_metrics_utils.apply_low_freq_fourier_filter
        res = tested_func(self.test_sig, 2)
        self.assertAlmostEqual(res[10].real, -0.0956296962962675, places=7)


class TestManney2011(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        # NOTE: this metric is a generator
        subset_data = self.data.sel(plev=slice(25000, 20000)).isel(
            time=slice(0, 1)
        )
        result = jetstream_metrics.manney_et_al_2011(subset_data)
        self.assertEqual(result["jet_core_id"].max(), 2)
        # jetstream_metrics.manney_et_al_2011(subset_data.transpose("lat", ...))


class TestPenaOrtiz2013(unittest.TestCase):
    def setUp(self):
        # TODO
        self.data = set_up_test_uv_data()

    def test_metric(self):
        # result = jetstream_metrics.penaortiz_et_al_2013(self.data)
        # self.assertIn("ws", result)
        pass

    def test_get_resultant_wind(self):
        pass

    def test_get_empty_local_wind_maxima(self):
        pass

    def test_num_timeunits_per_monthyear_with_local_wind_maxima(self):
        pass


class TestScreenSimmonds2013(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_metric(self):
        # result = jetstream_metrics.screen_and_simmonds_2013(self.data)
        pass


class TestKuang2014(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jetstream_metrics.kuang_et_al_2014
        lon_data = self.data.sel(plev=50000).isel(time=slice(0, 2))
        result = tested_func(lon_data)
        self.assertEqual(result["jet_ocurrence1_jet_centre2"].max(), 2)
        lon_data = self.data.sel(plev=slice(50000, 25000)).isel(
            time=slice(0, 2)
        )
        self.assertRaises(ValueError, lambda: tested_func(lon_data))


class TestFrancisVavrus2015(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        result = jetstream_metrics.francis_vavrus_2015(self.data)
        self.assertEqual(float(result["mci"].mean()), -0.01847001537680626)
        self.assertTrue(result["mci"].max() == 1)
        self.assertTrue(result["mci"].min() == -1)


class TestLocalWaveActivity(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_metric(self):
        result = jetstream_metrics.local_wave_activity(self.data)
        self.assertTrue(result)


class TestCattiaux2016(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_metric(self):
        # result = jetstream_metrics.cattiaux_et_al_2016(self.data)
        # self.assertTrue(result)
        pass


class TestBarnesSimpson2017(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        test_func = jetstream_metrics.barnes_simpson_2017
        result = test_func(self.data)
        self.assertEqual(float(result["north_atlantic_max_lats"].max()), 36.25)
        self.assertEqual(float(result["north_pacific_max_lats"].max()), 32.5)
        one_lon_data = self.data.isel(lon=slice(0, 1))
        only_north_atlantic_data = self.data.sel(lon=slice(280, 350))
        self.assertRaises(ValueError, lambda: test_func(one_lon_data))
        self.assertRaises(
            ValueError, lambda: test_func(only_north_atlantic_data)
        )


class TestGrisePolvani2017(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        result = jetstream_metrics.grise_polvani_2017(self.data)
        self.assertEqual(float(result["max_lat_0.01"].min()), 35.38)
        self.assertEqual(float(result["max_lat_0.01"].max()), 36.41)
        self.assertEqual(
            float(result["max_speed_0.01"].max()), 22.926441912353482
        )


class TestCeppi2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        result = jetstream_metrics.ceppi_et_al_2018(self.data)
        print(result)
        self.assertEqual(
            float(result["jet_lat_centroid"][0].data), 42.166801789276384
        )


# class TestKern2018(unittest.TestCase):
#     def setUp(self):
#         self.data = set_up_test_uv_data()

#     def test_metric(self):
#         result = jetstream_metrics.kern_et_al_2018(self.data)
#         self.assertTrue(result)


class TestBracegirdle2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        tested_func = jetstream_metrics.bracegirdle_et_al_2018
        test_data = self.data.sel(plev=85000)
        result = tested_func(test_data)
        self.assertRaises(ValueError, lambda: tested_func(self.data))
        self.assertEqual(float(result["seasonal_JPOS"].max()), 37.725)
        self.assertEqual(float(result["annual_JPOS"].max()), 37.725)
        self.assertEqual(round(float(result["seasonal_JSTR"].max()), 3), 8.589)
        self.assertEqual(round(float(result["annual_JSTR"].max()), 3), 8.589)


class TestJetStreamCoreIdentificationAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_ws_thresholds(self):
        tested_alg = (
            jetstream_metrics_utils.JetStreamCoreIdentificationAlgorithm
        )
        self.assertRaises(ValueError, lambda: tested_alg(self.data))
        test_data = self.data.isel(time=0, lon=0)
        self.assertRaises(ValueError, lambda: tested_alg(test_data, -10, 10))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, -10))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, 30))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, 10))

    def test_inner_funcs(self):
        tested_alg = (
            jetstream_metrics_utils.JetStreamCoreIdentificationAlgorithm
        )
        test_data = self.data.isel(time=0, lon=0)
        result = tested_alg(test_data)
        result.run()
        print(result._initial_core_ids.mean())
        self.assertEqual(result._initial_core_ids.mean(), 30.07608695652174)
        self.assertEqual(
            len(np.where(result._labelled_data["ws"] == "Core")[1]), 46
        )
        self.assertEqual(
            len(
                np.where(result._labelled_data["ws"] == "Potential Boundary")[
                    1
                ]
            ),
            81,
        )
        self.assertEqual(result.num_of_cores, 3)
        self.assertListEqual(
            result.final_jet_cores[0]["index_of_area"][0], [5, 15]
        )


class TestJetStreamOccurenceAndCentreAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_ws_thresholds(self):
        tested_alg = (
            jetstream_metrics_utils.JetStreamOccurenceAndCentreAlgorithm
        )
        test_data = self.data.isel(plev=0, time=0)
        self.assertRaises(ValueError, lambda: tested_alg(test_data, -10))

    def test_inner_functions(self):
        tested_alg = (
            jetstream_metrics_utils.JetStreamOccurenceAndCentreAlgorithm
        )
        test_data = self.data.isel(plev=4, time=0)
        result = tested_alg(test_data)
        result.run()
        self.assertEqual(
            float(result._jet_occurence["ws"].max()), 85.84358978271484
        )
        self.assertListEqual(result._jet_centres[0].tolist(), [5.0, 331.875])


if __name__ == "__main__":
    unittest.main()
