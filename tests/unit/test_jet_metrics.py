# -*- coding: utf-8 -*-

"""
    Tests for jet metrics (jetstream_metrics.py & jetstream_metrics_utils)
    Includes:
        – tests for each metric defined and expected outputs
        – test for each utility function used in the metric (e.g. core algorithms) 
"""

### imports
import unittest


### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"

# @pytest.mark.parametrize(["metric_name, all_metrics", ("'Woolings2010', ALL_METRICS", "this works")]) # example TODO
class TestKoch2006(unittest.TestCase):
    def test_basic(self):
        pass


class TestArcherCaldeira2008(unittest.TestCase):
    def test_basic(self):
        pass


class TestWoolings2010(unittest.TestCase):
    def test_basic(self):
        pass


class TestManney2011(unittest.TestCase):
    def test_basic(self):
        pass


class TestScreenSimmonds2013(unittest.TestCase):
    def test_basic(self):
        pass


class TestKuang2014(unittest.TestCase):
    def test_basic(self):
        pass


class TestFrancisVavrus2015(unittest.TestCase):
    def test_basic(self):
        pass


class TestLocalWaveActivity(unittest.TestCase):
    def test_basic(self):
        pass


class TestCattiaux2016(unittest.TestCase):
    def test_basic(self):
        pass


class TestCeppi2018(unittest.TestCase):
    def test_basic(self):
        pass


class TestKern2018(unittest.TestCase):
    def test_basic(self):
        pass


class TestSimpson2018(unittest.TestCase):
    def test_basic(self):
        pass


class TestChemkeMing2020(unittest.TestCase):
    def test_basic(self):
        pass


if __name__ == "__main__":
    unittest.main()