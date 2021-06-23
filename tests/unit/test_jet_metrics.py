# -*- coding: utf-8 -*-

"""
    Tests for jet metrics (jetstream_metrics.py & jetstream_metrics_utils)
    Includes:
        – tests for each metric defined and expected outputs
        – test for each utility function used in the metric (e.g. core algorithms) 
"""

### imports
import xarray as xr
import numpy as np
from metrics import jetstream_metrics, jetstream_metrics_utils, jetstream_metrics_dict
import unittest
from parameterized import parameterized


### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


MAX_VARIABLES = 4

def set_up_test_uv_data():
    u_data = xr.open_dataset("data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    v_data = xr.open_dataset("data/va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    data = xr.merge([u_data, v_data])
    data = data.sel(lat=slice(0, 90))
    data = data.isel(time=slice(0,100))
    return data


def set_up_test_u_data():
    data = xr.open_dataset("data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    data = data.sel(lat=slice(0, 90))
    data = data.isel(time=slice(0,100))
    return data


def set_up_test_zg_data():
    data = xr.open_dataset("data/zg_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    data = data.sel(lat=slice(0, 90))
    data = data.isel(time=slice(0,100))
    return data


class TestJetStreamMetricDict(unittest.TestCase): 
    def setUp(self):
        self.metric_dict = jetstream_metrics_dict.JETSTREAM_METRIC_DICT
        
    def test_metric_dict_keys(self):
        for metric_name in self.metric_dict.keys():
            self.assertIsInstance(metric_name, str)
    
    def test_metric_dict_values(self):
        for metric in self.metric_dict.values():
            self.assertIsInstance(metric, dict)
            self.assertEqual(len(metric.keys()),4)
            self.assertListEqual(list(metric.keys()), ["variables", "coords", "metric", "description"])

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

    @parameterized.expand([
    ("plev", 0, 100000),
    ("lat", -91, 91),
    ("lon", -1, 361),
   ])
    def test_each_coord(self, coord, min_value, max_value):
        for metric in self.metric_dict.values():
            if coord in metric["coords"].keys():
                self.assertEqual(len(metric["coords"][coord]), 2)
                self.assertGreaterEqual(min(metric["coords"][coord]), min_value)
                self.assertLessEqual(max(metric["coords"][coord]), max_value)

    def test_funcs(self):
        for metric in self.metric_dict.values():
            self.assertTrue(callable(metric['metric']))


class TestKoch2006(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()
    
    def test_basic(self):
        pass

    def test_get_all_plevs(self):
        ## make sure it returns an array
        self.assertIsInstance(jetstream_metrics_utils.get_all_plev_hPa(self.data), (np.ndarray))
        ## make sure it takes errors wrong types
        # self.assertRaises(TypeError, lambda: jetstream_metrics_utils.get_all_plev_hPa(['plev']))
        new_data = self.data.rename({'plev':'pl'})
        self.assertRaises(KeyError, lambda: jetstream_metrics_utils.get_all_plev_hPa(new_data))

    def test_sum_weighted_ws(self):    
        self.assertRaises(TypeError, lambda: jetstream_metrics_utils.get_sum_weighted_ws(self.data, 1))
        self.assertIsInstance(jetstream_metrics_utils.get_sum_weighted_ws(self.data, [0,100]), xr.DataArray)
        self.assertGreater(jetstream_metrics_utils.get_sum_weighted_ws(self.data, [0,1000]).max(), 0)

    def test_weighted_average_ws(self):
        self.assertRaises(TypeError, lambda: jetstream_metrics_utils.get_weighted_average_ws(self.data, 1))
        sum_weighted = jetstream_metrics_utils.get_sum_weighted_ws(self.data, [0,100])
        weighted_av = jetstream_metrics_utils.get_weighted_average_ws(sum_weighted, np.array([0,100]))
        self.assertGreater(weighted_av.min(), 0)
        self.assertGreaterEqual(weighted_av.max(), 0)

        weighted_av_threshold = weighted_av.where(weighted_av >= 30)
        self.assertGreater(weighted_av_threshold.max(), 0)
        self.assertGreaterEqual(weighted_av_threshold.min(), 0)



class TestArcherCaldeira2008(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


class TestWoolings2010(unittest.TestCase):
    def setUp(self):
        self.data  = set_up_test_u_data()
    
    def test_basic(self):
        pass


class TestManney2011(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


class TestScreenSimmonds2013(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_basic(self):
        pass


class TestKuang2014(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


class TestFrancisVavrus2015(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


class TestLocalWaveActivity(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_basic(self):
        pass


class TestCattiaux2016(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_basic(self):
        pass


class TestCeppi2018(unittest.TestCase):
    def setUp(self):
        self.data  = set_up_test_u_data()
    
    def test_basic(self):
        pass


class TestKern2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


class TestSimpson2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_basic(self):
        pass


class TestChemkeMing2020(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


class TestJetStreamOccurenceAndCentreAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


class TestJetStreamCoreIdentificationAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass

# @pytest.mark.parametrize(["metric_name, all_metrics", ("'Woolings2010', ALL_METRICS", "this works")]) # example TODO

if __name__ == "__main__":
    unittest.main()