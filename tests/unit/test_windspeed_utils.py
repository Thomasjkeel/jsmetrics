# -*- coding: utf-8 -*-

"""
    Tests for windspeed slice classes (windspeed_utils.py)
    Includes:
        – tests for all windspeed classes
"""

# imports
import unittest
from jsmetrics import windspeed_utils


class TestWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None, None)
        )
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None)
        )


class TestLatitudeWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None, None)
        )
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice("asf", None)
        )
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None, "dasf")
        )
        pass


class TestPressureLevelWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None, None)
        )
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice("asf", None)
        )
        self.assertRaises(
            TypeError, lambda: windspeed_utils.WindSpeedSlice(None, "dasf")
        )
        pass


if __name__ == "__main__":
    unittest.main()
