# -*- coding: utf-8 -*-

"""
    Tests for windspeed slice classes (windspeed_utils.py)
    Includes:
        – tests for all windspeed classes
"""

# imports
import unittest
from jsmetrics.utils import windspeed_utils
from . import set_up_test_uv_data


class TestGetZonalMean(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        tested_func = windspeed_utils.get_zonal_mean
        new_data = self.data.isel(time=0, plev=0)
        new_data = new_data.rename({"lon": "ln"})
        self.assertRaises(KeyError, lambda: tested_func(new_data))


class TestGetWindDirectionInDegrees(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        test_func = windspeed_utils.get_wind_direction_in_degrees
        self.assertEqual(round(test_func(u=10, v=20)), 207)


class TestWindSpeedSlice(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_basic(self):
        self.assertRaises(TypeError, lambda: windspeed_utils.WindSpeedSlice(None, None))
        self.assertRaises(TypeError, lambda: windspeed_utils.WindSpeedSlice(None))

    def test_get_values(self):
        wsslice = windspeed_utils.PressureLevelWindSpeedSlice(
            self.data.isel(time=0, plev=0)
        )
        self.assertEqual(wsslice.get_values().max(), 10)


class TestLatitudeWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(TypeError, lambda: windspeed_utils.WindSpeedSlice(None, None))
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice("asf", None)
        )
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None, "dasf")
        )


class TestPressureLevelWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(TypeError, lambda: windspeed_utils.WindSpeedSlice(None, None))
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice("asf", None)
        )
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None, "dasf")
        )


if __name__ == "__main__":
    unittest.main()
