# -*- coding: utf-8 -*-

"""
    Tests for general and data utils (data_utils.py)
"""

# imports
import cftime
from parameterized import parameterized
import unittest
import numpy as np
from jsmetrics.utils import data_utils
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


class TestAddNumDaysto360Datetime(unittest.TestCase):
    @parameterized.expand(
        [
            (
                cftime.Datetime360Day(day=1, month=1, year=2000, hour=1),
                360,
                cftime.Datetime360Day(day=1, month=1, year=2001, hour=1),
            ),
            (
                cftime.Datetime360Day(day=1, month=1, year=2000, hour=1),
                30,
                cftime.Datetime360Day(day=1, month=2, year=2000, hour=1),
            ),
            (
                cftime.Datetime360Day(day=1, month=1, year=2000, hour=1),
                29,
                cftime.Datetime360Day(day=30, month=1, year=2000, hour=1),
            ),
            (
                cftime.Datetime360Day(day=1, month=1, year=2000, hour=1),
                359,
                cftime.Datetime360Day(day=30, month=12, year=2000, hour=1),
            ),
        ]
    )
    def test_add_num_days_to360Datetime_expected_vals(
        self, test_date, days_to_add, expected_new_date
    ):
        added_test_date = data_utils.add_num_of_days_to_360Datetime(
            test_date, num_of_days_to_add=days_to_add
        )
        self.assertEqual(added_test_date, expected_new_date)

    def test_errors_raised(self):
        tested_func = data_utils.add_num_of_days_to_360Datetime
        test_date = cftime.Datetime360Day(day=1, month=1, year=2000, hour=1)
        self.assertRaises(
            AssertionError, lambda: tested_func(np.datetime64("1970-01-11"), 1)
        )
        self.assertRaises(
            AssertionError,
            lambda: tested_func(cftime.DatetimeAllLeap(day=1, month=1, year=2000), 1),
        )
        self.assertRaises(ValueError, lambda: tested_func(test_date, 0))
        self.assertRaises(ValueError, lambda: tested_func(test_date, -1))


class TestAddNumDaystoNoLeapDatetime(unittest.TestCase):
    @parameterized.expand(
        [
            (
                cftime.DatetimeNoLeap(day=1, month=1, year=2020, hour=1),
                365,
                cftime.DatetimeNoLeap(day=1, month=1, year=2021, hour=12),
            ),
            (
                cftime.DatetimeNoLeap(day=1, month=1, year=2019, hour=1),
                365,
                cftime.DatetimeNoLeap(day=1, month=1, year=2020, hour=12),
            ),
            (
                cftime.DatetimeNoLeap(day=28, month=2, year=2020, hour=1),
                1,
                cftime.DatetimeNoLeap(day=1, month=3, year=2020, hour=12),
            ),
        ]
    )
    def test_add_num_of_days_to_NoLeapDatetime(
        self, test_date, days_to_add, expected_new_date
    ):
        added_test_date = data_utils.add_num_of_days_to_NoLeapDatetime(
            test_date, num_of_days_to_add=days_to_add
        )
        self.assertEqual(added_test_date, expected_new_date)

    def test_errors_raised(self):
        tested_func = data_utils.add_num_of_days_to_NoLeapDatetime
        test_date = cftime.DatetimeNoLeap(day=1, month=1, year=2000, hour=1)
        self.assertRaises(
            AssertionError, lambda: tested_func(np.datetime64("1970-01-11"), 1)
        )
        self.assertRaises(
            AssertionError,
            lambda: tested_func(cftime.DatetimeAllLeap(day=1, month=1, year=2000), 1),
        )
        self.assertRaises(ValueError, lambda: tested_func(test_date, 0))
        self.assertRaises(ValueError, lambda: tested_func(test_date, -1))


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
