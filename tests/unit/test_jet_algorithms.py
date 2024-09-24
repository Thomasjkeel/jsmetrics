# -*- coding: utf-8 -*-

"""
    Tests for jet algorithms (jetstream_metrics.py & jetstream_metrics_utils)
    Includes:
        – tests for each algorithm defined and expected outputs
        – test for each utility function used in the metric
"""

# imports
import unittest
from jsmetrics import jet_core_algorithms

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
        self.assertEqual(round(float(result["jet_events_ws"].max()), 4), 8.7752)
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
        sub_data = self.data.isel(time=slice(0, 2), plev=slice(3, 6))
        result = tested_func(sub_data, ws_threshold=30)
        sub_data = self.data.isel(time=slice(0, 1), plev=slice(3, 6))
        result = tested_func(sub_data, ws_threshold=30)
        self.assertTrue(result)
        self.assertEqual(
            int(result["jet_occurence"].max("plev").sel(lat=50, lon=136.875)), 1
        )
        self.assertEqual(int(result["jet_occurence"].max("plev").sum()), 327)
        self.assertRaises(
            KeyError,
            lambda: tested_func(self.data.isel(time=0).drop_vars("time")),
        )


class TestManney2011(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.manney_et_al_2011
        subset_data = self.data.sel(plev=slice(50000, 10000)).isel(
            time=slice(0, 2), lon=slice(0, 100)
        )
        res = tested_func(subset_data, jet_core_plev_limit=(10000, 40000))
        self.assertEqual(int(res["jet_core_mask"].astype(int).max()), 1)
        subset_data = self.data.sel(plev=slice(50000, 10000)).isel(
            time=slice(0, 1), lon=slice(0, 100)
        )
        squeezed_data = subset_data.isel(plev=0)
        res = tested_func(squeezed_data, jet_core_plev_limit=(10000, 40000))
        res = tested_func(subset_data, jet_core_plev_limit=(10000, 40000))
        self.assertEqual(int(res["jet_core_mask"].astype(int).max()), 1)
        self.assertEqual(int(res["jet_region_mask"].astype(int).max()), 1)
        self.assertEqual(
            int(res["jet_core_mask"].isel(lon=0).where(lambda x: x).count()), 1
        )

        res = tested_func(
            subset_data, jet_core_plev_limit=(10000, 40000), jet_core_ws_threshold=10
        )
        self.assertEqual(
            int(res["jet_core_mask"].isel(lon=0).where(lambda x: x).count()), 5
        )
        self.assertRaises(
            KeyError,
            lambda: tested_func(
                self.data.isel(time=0).drop_vars("time"),
                jet_core_plev_limit=(10000, 40000),
            ),
        )
        self.assertRaises(
            KeyError,
            lambda: tested_func(self.data.isel(time=0), jet_core_plev_limit=False),
        )

    def test_check_diagonals(self):
        tested_func = jet_core_algorithms.manney_et_al_2011
        subset_data = self.data.sel(plev=slice(50000, 10000)).isel(
            time=slice(0, 1), lon=slice(0, 100)
        )
        res_diagonal = tested_func(
            subset_data, jet_core_plev_limit=(10000, 40000), check_diagonals=True
        )
        res_non_diagonal = tested_func(
            subset_data, jet_core_plev_limit=(10000, 40000), check_diagonals=False
        )

        self.assertNotEqual(
            int(res_non_diagonal["jet_core_mask"].where(lambda x: x).count()),
            int(res_diagonal["jet_core_mask"].where(lambda x: x).count()),
        )


class TestPenaOrtiz2013(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.penaortiz_et_al_2013
        subset_data = self.data.sel(plev=slice(25000, 20000)).isel(time=slice(0, 1))
        res = tested_func(subset_data)
        self.assertTrue("polar_front_jet" in res)
        self.assertEqual(res["subtropical_jet"].max(), 1)
        self.assertRaises(
            KeyError,
            lambda: tested_func(self.data.isel(time=0).drop_vars("time")),
        )

    def test_get_empty_local_wind_maxima(self):
        tested_func = jet_core_algorithms_components.get_empty_local_wind_maxima_data
        empty_local_wind_data = tested_func(self.data)
        self.assertEqual(empty_local_wind_data.max(), 0)
        self.assertTrue("local_wind_maxima" in empty_local_wind_data)

    def test_get_local_wind_maxima_by_timeunit(self):
        tested_func = jet_core_algorithms_components.get_local_wind_maxima_by_timeunit
        self.assertRaises(
            ValueError, lambda: tested_func(self.data.isel(time=0), ws_threshold=30)
        )


class TestKuang2014(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.kuang_et_al_2014
        lon_data = self.data.sel(plev=slice(50000, 25000)).isel(time=slice(0, 2))
        result = tested_func(lon_data)
        self.assertEqual(result["jet_centers"].max(), 1)
        self.assertEqual(result["jet_occurence"].max(), 1)
        plev_data = self.data.sel(plev=slice(25000, 25000)).isel(
            time=slice(0, 1), lat=slice(0, 10), lon=slice(0, 10)
        )
        tested_func(plev_data)
        self.assertRaises(
            KeyError,
            lambda: tested_func(self.data.isel(time=0).drop_vars("time")),
        )


class TestJetCoreIdentificationAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_core_algorithms.jet_core_identification_algorithm
        # NOTE: this metric is a generator
        subset_data = self.data.sel(plev=25000).isel(
            time=slice(0, 1), lat=slice(20, 30)
        )
        res = tested_func(subset_data)
        self.assertEqual(res["jet_core_id"].max(), 2)
        self.assertRaises(
            KeyError,
            lambda: tested_func(self.data.isel(time=0).drop_vars("time")),
        )


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
        test_data = self.data.isel(time=0, lon=0, lat=slice(20, 30))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, -10, 10))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, -10))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, 30))
        self.assertRaises(ValueError, lambda: tested_alg(test_data, 10, 10))

    def test_inner_funcs(self):
        tested_alg = jet_core_algorithms_components.JetStreamCoreIdentificationAlgorithm
        test_data = self.data.isel(time=0, lon=0, lat=slice(20, 50))
        result = tested_alg(test_data)
        result.run()
        print(result._initial_core_ids.mean())
        self.assertEqual(round(result._initial_core_ids.mean(), 3), 11.25)
        self.assertEqual(len(np.where(result._labelled_data["ws"] == "Core")[1]), 16)
        self.assertEqual(
            len(np.where(result._labelled_data["ws"] == "Potential Boundary")[1]),
            39,
        )
        self.assertEqual(result.num_of_cores, 1)
        self.assertListEqual(result.final_jet_cores[0]["index_of_area"][0], [8, 7])

    def test_classmethod(self):
        tested_alg = (
            jet_core_algorithms_components.JetStreamCoreIdentificationAlgorithm.run_algorithm
        )
        test_data = self.data.isel(time=0, lon=0, lat=slice(20, 30))
        result = tested_alg(test_data)
        self.assertEqual(round(result._initial_core_ids.mean(), 3), 7.75)
        self.assertEqual(len(np.where(result._labelled_data["ws"] == "Core")[1]), 2)


class TestJetStreamOccurenceAndCentreAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data().isel(time=slice(0, 1))
        self.tested_alg = (
            jet_core_algorithms_components.JetStreamOccurenceAndCentreAlgorithm
        )

    def test_ws_thresholds(self):
        test_data = self.data.isel(plev=0, time=0, lat=slice(20, 30), lon=slice(0, 10))
        self.assertRaises(ValueError, lambda: self.tested_alg(test_data, -10))

    def test_inner_functions(self):
        test_data = self.data.isel(plev=4, time=0, lat=slice(20, 30), lon=slice(0, 10))
        result = self.tested_alg(test_data)
        result.run()
        self.assertEqual(round(float(result._jet_occurence["ws"].max()), 3), 69.214)
        self.assertListEqual(result._jet_centres[0].tolist(), [30.0, 16.875])

    def test_cls_method(self):
        test_data = self.data.isel(plev=4, lat=slice(20, 30), lon=slice(0, 10), time=0)
        result = self.tested_alg.run_algorithm(test_data)
        self.assertEqual(round(float(result._jet_occurence["ws"].max()), 3), 69.214)


if __name__ == "__main__":
    unittest.main()
