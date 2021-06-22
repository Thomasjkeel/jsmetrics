# -*- coding: utf-8 -*-

"""
    Tests for windspeed slice classes (windspeed_utils.py)
    Includes:
        – tests for all windspeed classes
"""

### imports
import numpy as np
import unittest
from metrics import windspeed_utils
# from parameterized import parameterized
# import math

class TestWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice(None, None))
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice('asf', None))
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice(None, 'dasf'))

#     @parameterized.expand([
#     ("negative", -1.5, -2.0),
#     ("integer", 1, 1.0),
#     ("large fraction", 1.6, 1),
#    ])
#     def test_floor(self, name, input, expected):
#         self.assertEqual(math.floor(input), expected)


class TestLatitudeWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice(None, None))
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice('asf', None))
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice(None, 'dasf'))


class TestPressureLevelWindSpeedSlice(unittest.TestCase):
    def test_basic(self):
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice(None, None))
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice('asf', None))
        self.assertRaises(ValueError, lambda: windspeed_utils.WindSpeedSlice(None, 'dasf'))


if __name__ == "__main__":
    unittest.main()