# -*- coding: utf-8 -*-

"""
    Tests for computing metrics (compute_metrics.py)
    Includes:
        – methods for subsetting
        – checks of data before computing metric 
"""

### imports
from tests import unit
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
        UKESM1_SSP585_U = xr.open_dataset("data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
        UKESM1_SSP585_V = xr.open_dataset("data/va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
        UKESM1_SSP585 = xr.merge([UKESM1_SSP585_U, UKESM1_SSP585_V])
        ukesm1_ssp585 = compute_metrics.MetricComputer(UKESM1_SSP585)
        self.data = ukesm1_ssp585.sel(lat=slice(0, 90))
        self.data = self.data.isel(time=slice(0,100))
        self.all_metrics = JETSTREAM_METRIC_DICT

    def test_metric_class_inputs(self):
        # compute_metrics.MetricComputer()
        pass

    def test_available_metrics(self):
        pass

    def test_get_variables(self):
        pass

    def test_swap_coords(self):
        pass
 
    def test_sel(self):
        with self.assertRaises(TypeError): #TODO
            self.data.sel(fake=slice(0, 90))
            
    def test_subset_data(self):
        pass

    def test_all_metrics(self):
        pass

    def test_check_all_variables(self):
        pass

    def test_check_all_chords(self):
        pass

    def test_check_coords_meet_reqs(self):
        pass

    def test_compute_metrics(self):
        pass


if __name__ == "__main__":
    unittest.main()