# -*- coding: utf-8 -*-

"""
    Tests for computing metrics (compute_metrics.py)
    Includes:
        – methods for subsetting
        – checks of data before computing metric 
"""

### imports
import numpy as np
import xarray as xr
from metrics import compute_metrics
from . import set_up_test_uv_data, set_up_test_u_data, set_up_test_zg_data, set_up_nan_dataset
from metrics.jetstream_metrics_dict import JETSTREAM_METRIC_DICT
import unittest

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"



class TestMetricComputer(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()
        
    def test_inputs(self):
        theClass = compute_metrics.MetricComputer
        self.assertRaises(AssertionError, lambda: theClass("wrong", {"a":1}))
        self.assertRaises(AssertionError, lambda: theClass(self.data, {}))
        self.assertRaises(AssertionError, lambda: theClass(self.data, [0]))
        self.assertRaises(AssertionError, lambda: theClass(self.data['ua'], {"a":1}))
    
    def test_inputs_with_available_metrics(self):
        theClassMethod = compute_metrics.MetricComputer.with_available_metrics
        metric_computer = theClassMethod(self.data, all_metrics=JETSTREAM_METRIC_DICT)
        self.assertTrue(hasattr(metric_computer, 'available_metrics'))
        self.assertIsInstance(metric_computer.available_metrics, list)
        self.assertRaises(TypeError, lambda: theClassMethod(self.data, {"a":1}))

    def test_basic(self):
        metric_computer = compute_metrics.MetricComputer(self.data, all_metrics=JETSTREAM_METRIC_DICT)
        self.assertTrue(hasattr(metric_computer, 'data'))
        self.assertTrue(hasattr(metric_computer, 'variable_list'))
        self.assertTrue(hasattr(metric_computer, 'all_metrics'))
        self.assertIsInstance(metric_computer.variable_list, list)
        self.assertIsInstance(metric_computer.all_metrics, dict)
    
    def test_available_metrics(self):
        metric_computer = compute_metrics.MetricComputer(self.data, all_metrics=JETSTREAM_METRIC_DICT)
        metric_computer.get_available_metrics()
        self.assertTrue(hasattr(metric_computer, 'available_metrics'))
        for metric in metric_computer.available_metrics:
            self.assertTrue(list(metric.keys())[0] in JETSTREAM_METRIC_DICT.keys())
            self.assertTrue(len(metric) == 1)
            self.assertIsInstance(list(metric.values())[0], str)

    def test_get_variables(self):
        metric_computer = compute_metrics.MetricComputer(self.data, all_metrics=JETSTREAM_METRIC_DICT)
        self.assertListEqual(metric_computer.variable_list, ['ua', 'va'])

    def test_swap_coords(self):
        pass
 
    def test_sel(self):
        metric_computer = compute_metrics.MetricComputer(self.data, all_metrics=JETSTREAM_METRIC_DICT)
        self.assertRaises(ValueError, lambda: metric_computer.sel(fake=slice(0, 90)))
        new1 = metric_computer.sel(lon=slice(0, 90))
        new2 = self.data.sel(lon=slice(0, 90))
        self.assertListEqual(list(new1.data.lon.values), list(new2.lon.values))

    def test_isel(self):
        metric_computer = compute_metrics.MetricComputer(self.data, all_metrics=JETSTREAM_METRIC_DICT)
        self.assertRaises(ValueError, lambda: metric_computer.isel(fake=slice(0, 90)))
        new1 = metric_computer.isel(lon=slice(0, 90))
        new2 = self.data.isel(lon=slice(0, 90))
        self.assertListEqual(list(new1.data.lon.values), list(new2.lon.values))

    def test_compute_metric_from_data(self):
        metric_computer = compute_metrics.MetricComputer.with_available_metrics(self.data, all_metrics=JETSTREAM_METRIC_DICT)
        test_metric_name = 'FrancisVavrus2015'
        result = metric_computer.compute_metric_from_data(test_metric_name)
        self.assertEqual(float(result['mci'][0][0][0]), -0.015743955969810486)
        bad_metric_computer = compute_metrics.MetricComputer(self.data, all_metrics={"1":1})
        self.assertRaises(KeyError, lambda: bad_metric_computer.compute_metric_from_data(test_metric_name))
        #TODO: add subset and calc kwarg tests
 

class TestComputeMetricFunctions(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_subset_data(self):
        tested_func = compute_metrics.subset_data
        self.assertRaises(AssertionError, lambda: tested_func([], 'hello'))
        test_metric = JETSTREAM_METRIC_DICT['FrancisVavrus2015']
        result = tested_func(self.data, test_metric)
        plev_coords_before = list(self.data['plev'].values)
        plev_coords_after = [result['plev'].values.tolist()]
        self.assertTrue(plev_coords_after[0] in plev_coords_before)
        self.assertNotEqual(len(plev_coords_before), len(plev_coords_after))

    def test_subset_data_with_ignore_coords(self):
        tested_func = compute_metrics.subset_data
        test_metric = JETSTREAM_METRIC_DICT['FrancisVavrus2015']
        self.assertRaises(AssertionError, lambda: tested_func(self.data, test_metric, ignore_coords="wrong"))
        self.assertRaises(AssertionError, lambda: tested_func(self.data, test_metric, ignore_coords=("wrong")))
        result = tested_func(self.data, test_metric, ignore_coords=['plev'])
        plev_coords_before = list(self.data['plev'].values)
        plev_coords_after = result['plev'].values.tolist()
        self.assertEqual(len(plev_coords_before), len(plev_coords_after))

    def test_get_coords_to_subset(self):
        tested_func = compute_metrics.get_coords_to_subset
        pass

    def test_check_all_variables(self):
        tested_func = compute_metrics.check_all_variables_available
        pass

    def test_check_all_chords(self):
        tested_func = compute_metrics.check_all_coords_available


    def test_check_coords_meet_reqs(self):
        tested_func = compute_metrics.check_if_coord_vals_meet_reqs

    def test_compute_metrics(self):
        tested_func = compute_metrics.compute_metric


    def test_swap_coords(self):
        tested_func = compute_metrics.swap_coord_order


if __name__ == "__main__":
    unittest.main()