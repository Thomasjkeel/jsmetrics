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
from metrics.jetstream_metrics_dict import JETSTREAM_METRIC_DICT
import unittest

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


class TestMetricComputer(unittest.TestCase):
    """
        maybe split these classes
    """
    def setUp(self):
        u_data = xr.open_dataset("data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
        v_data = xr.open_dataset("data/va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
        uv_data = xr.merge([u_data, v_data])
        self.data = compute_metrics.MetricComputer(uv_data)
        self.data = self.data.sel(lat=slice(0, 90))
        self.data = self.data.isel(time=slice(0,100))
        self.all_metrics = JETSTREAM_METRIC_DICT

    def test_basic(self):
        self.assertRaises(ValueError, lambda: compute_metrics.MetricComputer(None))
        self.assertRaises(ValueError, lambda: compute_metrics.MetricComputer('asf'))
        self.assertRaises(ValueError, lambda: compute_metrics.MetricComputer((None, 'dasf')))

    def test_available_metrics(self):
        pass

    def test_get_variables(self):
        pass

    def test_swap_coords(self):
        pass
 
    def test_sel(self):
        with self.assertRaises(ValueError): #TODO
            self.data.sel(fake=slice(0, 90))


class TestComputeMetricFunctions(unittest.TestCase):
    def test_subset_data(self):
        pass

    def test_all_metrics(self):
        """
            maybe not needed as in test_jet_metrics.py TODO: remove if not needed
        """
        pass

    def test_check_all_variables(self):
        pass

    def test_check_all_chords(self):
        pass

    def test_check_coords_meet_reqs(self):
        pass

    def test_compute_metrics(self):
        pass

    def test_swap_coords(self):
        pass


if __name__ == "__main__":
    unittest.main()