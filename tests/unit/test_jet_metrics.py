# -*- coding: utf-8 -*-

"""
    Tests for jet metrics (jetstream_metrics.py & jetstream_metrics_utils)
    Includes:
        – tests for each metric defined and expected outputs
        – test for each utility function used in the metric (e.g. core algorithms) 
"""

### imports
import xarray as xr
from metrics import jetstream_metrics
import unittest


### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


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


class TestKoch2006(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        pass


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