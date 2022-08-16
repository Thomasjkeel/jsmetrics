# -*- coding: utf-8 -*-

"""
    Tests for spatial utils (spatial_utils.py)
"""

# imports
# import shapely.geometry

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


import unittest
import numpy as np
from jsmetrics import spatial_utils
from . import set_up_test_uv_data, set_up_test_zg_data


class TestStandardiseDiffsByMakingAllMostCommonDiff(unittest.TestCase):
    def test_basic(self):
        test_func = spatial_utils._standardise_diffs_by_making_all_most_common_diff
        test_func([0, 1])


class TestQuadrant_area(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(
            ValueError,
            lambda: spatial_utils._quadrant_area(np.array([0, 1, 2]), None, None),
        )


class TestCalcSpatialIntegral(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        test_func = spatial_utils.calc_spatial_integral
        res = test_func(
            self.data.isel(time=0, plev=0)["ua"],
            lon_name="lon",
            lat_name="lat",
        )
        self.assertEqual(round(float(res)), -358693799460118)


class TestCalcSpatialMean(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        test_func = spatial_utils.calc_spatial_mean
        res = test_func(
            self.data.isel(time=0, plev=0)["ua"],
            lon_name="lon",
            lat_name="lat",
        )
        self.assertEqual(round(float(res)), -1)


class TestCalcTotalGreatCircleDistanceAlongLine(unittest.TestCase):
    def test_basic(self):
        test_func = spatial_utils.calc_total_great_circle_distance_along_line
        self.assertTrue(np.isnan(test_func([0, 1])))


class TestGetOneContourLinestring(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_basic(self):
        test_func = spatial_utils.get_one_contour_linestring
        test_data = self.data.isel(time=0, plev=0)["zg"]
        test_func(test_data.isel(lat=slice(0, 6), lon=slice(0, 5)), 100)
        # leaves the work-around for valueerrors untested (because requires lots of data to find edge case)
