# -*- coding: utf-8 -*-

"""
    Tests for general and data utils (data_utils.py)
"""

# imports
import unittest

# from parameterized import parameterized
from jsmetrics import data_utils
from . import (
    set_up_test_uv_data,
)

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


class TestXarrayDataCheck(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_check_coords_in_data(self):
        tested_func = data_utils.check_coords_in_data
        self.assertRaises(KeyError, lambda: tested_func(self.data, ("notthere",)))

    def test_check_variables_in_data(self):
        tested_func = data_utils.check_variables_in_data
        self.assertRaises(
            KeyError,
            lambda: tested_func(self.data, req_variables=["notthere"]),
        )

    def test_check_at_least_n_plevs_in_data(self):
        tested_func = data_utils.check_at_least_n_plevs_in_data
        tested_func(self.data, n_plevs=1)
        self.assertRaises(ValueError, lambda: tested_func(self.data, n_plevs=10))

    def test_check_if_data_is_xarray_datatype(self):
        tested_func = data_utils.check_if_data_is_xarray_datatype
        tested_func(self.data)
        new_data = self.data["ua"].data
        self.assertRaises(TypeError, lambda: tested_func(new_data))


class TestLocalMinimaMaxima(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_get_local_maxima(self):
        tested_func = data_utils.get_local_maxima
        self.assertListEqual(
            list(tested_func(self.data["ua"].isel(time=0, plev=0, lon=0).data)[0]),
            [4, 30, 38, 44, 56, 61, 70],
        )

    def test_get_local_minima(self):
        tested_func = data_utils.get_local_minima
        self.assertListEqual(
            list(tested_func(self.data["ua"].isel(time=0, plev=0, lon=0).data)[0]),
            [43, 48, 58, 65],
        )


class TestGetNumOfDecimalPlaces(unittest.TestCase):
    def test_func(self):
        tested_func = data_utils.get_num_of_decimal_places
        self.assertEqual(tested_func(2.33), 2)
        self.assertEqual(tested_func(2), 0)
