# -*- coding: utf-8 -*-

"""
    Tests for general utils (data_utils.py)
    Includes:
    TODO
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
        self.assertRaises(
            KeyError, lambda: tested_func(self.data, ("notthere",))
        )

    def test_check_variables_in_data(self):
        tested_func = data_utils.check_variables_in_data
        self.assertRaises(
            KeyError,
            lambda: tested_func(self.data, req_variables=["notthere"]),
        )

    def test_check_at_least_n_plevs_in_data(self):
        tested_func = data_utils.check_at_least_n_plevs_in_data
        self.assertRaises(ValueError, lambda: tested_func(self.data))

    def test_check_if_data_is_xarray_datatype(self):
        tested_func = data_utils.check_if_data_is_xarray_datatype
        tested_func(self.data)
        new_data = self.data["ua"].data
        self.assertRaises(TypeError, lambda: tested_func(new_data))
