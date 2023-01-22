# -*- coding: utf-8 -*-

"""
    Tests for jet algorithms (jetstream_metrics.py & jetstream_metrics_utils)
    Includes:
        – tests for each algorithm defined and expected outputs
        – test for each utility function used in the metric
"""

# imports
import unittest
from jsmetrics.metrics import jet_core_algorithms

# from parameterized import parameterized
import xarray as xr
import numpy as np
from jsmetrics.metrics import (
    jet_core_algorithms_components,
)
from . import set_up_test_uv_data

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


class TestKoch2006(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.koch_et_al_2006
        result = tested_func(self.data, ws_threshold=8)
        # check an exact value. Is this necessary?
        self.assertIsInstance(result, xr.Dataset)
        self.assertEqual(float(result["weighted_average_ws"].max()), 8.775158882141113)
        new_data = self.data.isel(plev=0)
        self.assertRaises(ValueError, lambda: tested_func(new_data))

    def test_get_all_hPa_list(self):
        tested_func = jet_core_algorithms_components.get_all_hPa_list
        # make sure it returns an array
        self.assertIsInstance(tested_func(self.data), (np.ndarray))
        # make sure it takes errors wrong types
        # self.assertRaises(TypeError,\
        # lambda: jetstream_metrics_utils.get_all_plev_hPa(['plev']))
        new_data = self.data.rename({"plev": "pl"})
        self.assertRaises(KeyError, lambda: tested_func(new_data))
        new_data2 = self.data.copy()
        new_data2["plev"] = self.data.plev.assign_attrs(units="wrong")
        self.assertRaises(ValueError, lambda: tested_func(new_data2))

    def test_sum_weighted_ws(self):
        tested_func = jet_core_algorithms_components.get_sum_weighted_ws
        self.assertRaises(TypeError, lambda: tested_func(self.data, 1))
        new_data = self.data.rename({"plev": "pl"})
        self.assertRaises(KeyError, lambda: tested_func(new_data, [10, 20]))
        sum_weighted = tested_func(self.data, [0, 100])
        self.assertIsInstance(sum_weighted, xr.DataArray)
        self.assertGreater(sum_weighted.max(), 0)

    def test_weighted_average_ws(self):
        tested_func = jet_core_algorithms_components.get_weighted_average_ws

        self.assertRaises(TypeError, lambda: tested_func(self.data, 1))
        sum_weighted = jet_core_algorithms_components.get_sum_weighted_ws(
            self.data, [0, 100]
        )
        weighted_av = tested_func(sum_weighted, np.array([0, 100]))
        self.assertGreater(weighted_av.min(), 0)
        self.assertGreaterEqual(weighted_av.max(), 0)

        weighted_av_threshold = weighted_av.where(weighted_av >= 30)
        self.assertGreater(weighted_av_threshold.max(), 0)
        self.assertGreaterEqual(weighted_av_threshold.min(), 0)


class TestSchiemann2009(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.schiemann_et_al_2009
        sub_data = self.data.isel(time=slice(0, 1), lon=slice(0, 30))
        result = tested_func(sub_data)
        self.assertTrue(result)


class TestManney2011(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        test_func = jet_core_algorithms.manney_et_al_2011
        # NOTE: this metric is a generator
        subset_data = self.data.sel(plev=slice(25000, 20000)).isel(time=slice(0, 1))
        res = test_func(subset_data.sel(plev=25000))
        self.assertEqual(res["jet_core_id"].max(), 2)


class TestPenaOrtiz2013(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.penaortiz_et_al_2013
        subset_data = self.data.sel(plev=slice(25000, 20000)).isel(time=slice(0, 1))
        res = tested_func(subset_data)
        self.assertTrue("polar_front_jet" in res)
        self.assertEqual(res["subtropical_jet"].max(), 1)

    def test_get_empty_local_wind_maxima(self):
        test_func = jet_core_algorithms_components.get_empty_local_wind_maxima_data
        empty_local_wind_data = test_func(self.data)
        self.assertEqual(empty_local_wind_data.max(), 0)
        self.assertTrue("local_wind_maxima" in empty_local_wind_data)

    def test_get_local_wind_maxima_by_timeunit(self):
        test_func = jet_core_algorithms_components.get_local_wind_maxima_by_timeunit
        self.assertRaises(ValueError, lambda: test_func(self.data.isel(time=0)))


class TestKuang2014(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.kuang_et_al_2014
        lon_data = self.data.sel(plev=50000).isel(time=slice(0, 2))
        result = tested_func(lon_data)
        self.assertEqual(result["jet_ocurrence1_jet_centre2"].max(), 2)
        lon_data = self.data.sel(plev=slice(50000, 25000)).isel(time=slice(0, 2))
        self.assertRaises(ValueError, lambda: tested_func(lon_data))
        plev_data = self.data.sel(plev=slice(25000, 25000)).isel(time=slice(0, 2))
        tested_func(plev_data)


class TestJetStreamCoreIdentificationAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_repr(self):
        tested_alg = jet_core_algorithms_components.JetStreamCoreIdentificationAlgorithm
        test_data = self.data.isel(time=0, lon=0)
        jetCoreAlg = tested_alg(test_data)

        repr(jetCoreAlg)
        jetCoreAlg.run()
        repr(jetCoreAlg)

    def test_ws_thresholds(self):
        tested_alg = jet_core_algorithms_components.JetStreamCoreIdentificationAlgorithm
        self.assertRaises(ValueError, lambda: tested_alg(self.data))
        test_data = self.data.isel(time=0, lon=0)
        self.assertRaises(ValueError, lambda: tested_alg(test_data, -10, 10))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, -10))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, 30))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, 10))

    def test_inner_funcs(self):
        tested_alg = jet_core_algorithms_components.JetStreamCoreIdentificationAlgorithm
        test_data = self.data.isel(time=0, lon=0)
        result = tested_alg(test_data)
        result.run()
        print(result._initial_core_ids.mean())
        self.assertEqual(result._initial_core_ids.mean(), 30.07608695652174)
        self.assertEqual(len(np.where(result._labelled_data["ws"] == "Core")[1]), 46)
        self.assertEqual(
            len(np.where(result._labelled_data["ws"] == "Potential Boundary")[1]),
            81,
        )
        self.assertEqual(result.num_of_cores, 3)
        self.assertListEqual(result.final_jet_cores[0]["index_of_area"][0], [16, 4])

    def test_classmethod(self):
        tested_alg = (
            jet_core_algorithms_components.JetStreamCoreIdentificationAlgorithm.run_algorithm
        )
        test_data = self.data.isel(time=0, lon=0)
        result = tested_alg(test_data)
        self.assertEqual(result._initial_core_ids.mean(), 30.07608695652174)
        self.assertEqual(len(np.where(result._labelled_data["ws"] == "Core")[1]), 46)


class TestJetStreamOccurenceAndCentreAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data().isel(time=slice(0, 1))
        self.tested_alg = (
            jet_core_algorithms_components.JetStreamOccurenceAndCentreAlgorithm
        )

    def test_ws_thresholds(self):
        test_data = self.data.isel(plev=0, time=0)
        self.assertRaises(ValueError, lambda: self.tested_alg(test_data, -10))

    def test_inner_functions(self):
        test_data = self.data.isel(plev=4, time=0)
        result = self.tested_alg(test_data)
        result.run()
        self.assertEqual(float(result._jet_occurence["ws"].max()), 85.84358978271484)
        self.assertListEqual(result._jet_centres[0].tolist(), [0.0, 247.5])

    def test_cls_method(self):
        test_data = self.data.isel(plev=4, time=0)
        result = self.tested_alg.run_algorithm(test_data)
        self.assertEqual(float(result._jet_occurence["ws"].max()), 85.84358978271484)


if __name__ == "__main__":
    unittest.main()
