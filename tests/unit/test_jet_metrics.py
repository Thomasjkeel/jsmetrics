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
from . import set_up_test_uv_data, set_up_test_u_data, set_up_test_zg_data, set_up_nan_dataset, make_fake_seasonal_data, make_fake_data
import unittest
from parameterized import parameterized


### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


MAX_VARIABLES = 4

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
    
    def test_metric(self):
        result = jetstream_metrics.koch_et_al_2006(self.data, ws_threshold=8)
        ## check an exact value. Is this necessary?
        self.assertIsInstance(result, xr.Dataset)
        self.assertEqual(float(result['weighted_average_ws'].max()), 8.775158882141113)

    def test_get_all_plevs(self):
        tested_func = jetstream_metrics_utils.get_all_plev_hPa
        ## make sure it returns an array
        self.assertIsInstance(tested_func(self.data), (np.ndarray))
        ## make sure it takes errors wrong types
        # self.assertRaises(TypeError, lambda: jetstream_metrics_utils.get_all_plev_hPa(['plev']))
        new_data = self.data.rename({'plev':'pl'})
        self.assertRaises(KeyError, lambda: tested_func(new_data))

    def test_sum_weighted_ws(self):    
        tested_func = jetstream_metrics_utils.get_sum_weighted_ws
        self.assertRaises(TypeError, lambda: tested_func(self.data, 1))
        new_data = self.data.rename({'plev':'pl'})
        self.assertRaises(KeyError, lambda: tested_func(new_data, [10,20]))
        sum_weighted = tested_func(self.data, [0,100])
        self.assertIsInstance(sum_weighted, xr.DataArray)
        self.assertGreater(sum_weighted.max(), 0)

    def test_weighted_average_ws(self):
        tested_func = jetstream_metrics_utils.get_weighted_average_ws

        self.assertRaises(TypeError, lambda: tested_func(self.data, 1))
        sum_weighted = jetstream_metrics_utils.get_sum_weighted_ws(self.data, [0,100])
        weighted_av = tested_func(sum_weighted, np.array([0,100]))
        self.assertGreater(weighted_av.min(), 0)
        self.assertGreaterEqual(weighted_av.max(), 0)

        weighted_av_threshold = weighted_av.where(weighted_av >= 30)
        self.assertGreater(weighted_av_threshold.max(), 0)
        self.assertGreaterEqual(weighted_av_threshold.min(), 0)


class TestArcherCaldeira2008(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        # result = jetstream_metrics.archer_caldeira_2008(self.data)
        pass


class TestWoolings2010(unittest.TestCase):
    def setUp(self):
        self.data  = set_up_test_u_data()
        self.data = make_fake_seasonal_data(self.data)
    
    def test_metric(self):
        result = jetstream_metrics.woolings_et_al_2010(self.data, filter_freq=1, window_size=2)
        self.assertIsInstance(result, xr.Dataset)
        print(result['filtered_max_ws'][0])
        self.assertEqual(result['filtered_max_lats'][0], 36.25)
        self.assertEqual(result['filtered_max_ws'][0], 44.532657623291016)

    def test_get_zonal_mean(self):
        tested_func = jetstream_metrics_utils.get_zonal_mean
        new_data = self.data.rename({'lon':'ln'})
        self.assertRaises(KeyError, lambda: tested_func(new_data))
        self.assertIsInstance(tested_func(self.data), xr.Dataset)

    def test_apply_lancoz_filter(self):
        tested_func = jetstream_metrics_utils.apply_lancoz_filter
        self.assertRaises(AssertionError, lambda: tested_func(self.data, -2, 1))
        self.assertRaises(AssertionError, lambda: tested_func(self.data, 2, -1))
        self.assertRaises(AssertionError, lambda: tested_func(self.data, 2, 1))
        self.assertRaises(AssertionError, lambda: tested_func(self.data, self.data['time'].count()+2, 1))
        self.assertRaises(AssertionError, lambda: tested_func(self.data, 2, self.data['time'].count()+1))
        self.assertEqual(float(tested_func(self.data, 2,4).max()), 99.514892578125)

    def test_get_latitude_and_speed_where_max_ws(self):
        tested_func = jetstream_metrics_utils.get_latitude_and_speed_where_max_ws
        self.assertRaises(AttributeError, lambda: tested_func(['lol']))
        tested_data = self.data['ua'].isel(plev=0, lon=0, time=0)
        self.assertEqual(tested_func(tested_data)[0], 70.0)
        self.assertEqual(tested_func(tested_data)[1], 3.105090856552124)
        self.assertRaises(KeyError, lambda: tested_func(tested_data.rename({'lat':'lt'})))
        nan_dataset = set_up_nan_dataset()
        self.assertEqual(tested_func(nan_dataset), (None, None))

    def test_apply_fourier_filter(self):
        tested_func = jetstream_metrics_utils.apply_low_freq_fourier_filter
                

class TestManney2011(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        ## NOTE: this metric is a generator
        result = jetstream_metrics.manney_et_al_2011(self.data)
        self.assertRaises(ValueError, lambda: next(result))
        lon_data = self.data.isel(lon=0)
        result = jetstream_metrics.manney_et_al_2011(lon_data)
        current = next(result)
        self.assertEqual(current.core_ids.mean(), 30.07608695652174)
        self.assertEqual(len(np.where(current.output['ws'] == 'Core')[1]), 46)
        self.assertEqual(len(np.where(current.output['ws'] == 'Potential Boundary')[1]), 81)
        alg_results = current.run()
        self.assertEqual(len(alg_results), 3)
        self.assertListEqual(alg_results[0]['index_of_area'][0], [5,15])


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
        result = jetstream_metrics.kuang_et_al_2014(self.data)
        self.assertRaises(ValueError, lambda: next(result))
        lon_data = self.data.sel(plev=50000)
        result = jetstream_metrics.kuang_et_al_2014(lon_data)
        current = next(result)
        self.assertEqual(float(current.jet_occurence['ws'].max()), 50.63158416748047) 
        current.run()
        self.assertListEqual(current.jet_centres[0].tolist(), [28.75, 60.])


class TestFrancisVavrus2015(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        result = jetstream_metrics.francis_vavrus_2015(self.data)
        self.assertEqual(float(result['mci'].mean()), -0.01847001537680626)
        self.assertTrue(result['mci'].max() == 1)
        self.assertTrue(result['mci'].min() == -1)


class TestLocalWaveActivity(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_metric(self):
        # result = jetstream_metrics.local_wave_activity(self.data)
        pass


class TestCattiaux2016(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()
        
    def test_metric(self):
        # result = jetstream_metrics.cattiaux_et_al_2016(self.data)
        pass


class TestCeppi2018(unittest.TestCase):
    def setUp(self):
        self.data  = set_up_test_u_data()
    
    def test_metric(self):
        result = jetstream_metrics.ceppi_et_al_2018(self.data)
        current = next(result)
        self.assertEqual(current, 37.96227342516621)

class TestKern2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        # result = jetstream_metrics.kern_et_al_2018(self.data)
        pass

class TestSimpson2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_metric(self):
        # result = jetstream_metrics.simpson_et_al_2018(self.data)
        pass


class TestChemkeMing2020(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        # result = jetstream_metrics.chemke_and_ming_2020(self.data)
        pass


class TestJetStreamCoreIdentificationAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_ws_thresholds(self):
        tested_alg = jetstream_metrics_utils.JetStreamCoreIdentificationAlgorithm
        self.assertRaises(ValueError, lambda: tested_alg(self.data, 40, 30))
        test_data = self.data.isel(lon=0, time=0)
        self.assertRaises(AssertionError, lambda: tested_alg(test_data,-10,10))
        self.assertRaises(AssertionError, lambda: tested_alg(test_data,10,-10))
        self.assertRaises(AssertionError, lambda: tested_alg(test_data,10,30))
        self.assertRaises(AssertionError, lambda: tested_alg(test_data,10,10))



class TestJetStreamOccurenceAndCentreAlgorithm(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_ws_thresholds(self):
        tested_alg = jetstream_metrics_utils.JetStreamOccurenceAndCentreAlgorithm
        self.assertRaises(ValueError, lambda: tested_alg(self.data, 40))
        test_data = self.data.isel(plev=0, time=0)
        self.assertRaises(AssertionError, lambda: tested_alg(test_data,-10))



if __name__ == "__main__":
    unittest.main()