# imports
import unittest
from jsmetrics.metrics import waviness_metrics
from . import set_up_test_uv_data, set_up_test_zg_data

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


class TestFrancisVavrus2015(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        result = waviness_metrics.francis_vavrus_2015(self.data)
        self.assertEqual(float(result["mci"].mean()), -0.01847001537680626)
        self.assertTrue(result["mci"].max() == 1)
        self.assertTrue(result["mci"].min() == -1)


# class TestLocalWaveActivity(unittest.TestCase):
#     def setUp(self):
#         self.data = set_up_test_zg_data()

#     def test_metric(self):
#         result = waviness_metrics.local_wave_activity(self.data)
#         self.assertTrue(result)


class TestCattiaux2016(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_metric(self):
        tested_func = waviness_metrics.cattiaux_et_al_2016
        subset_data = self.data.isel(plev=0)
        res = tested_func(subset_data)
        tested_func(subset_data["zg"])
        self.assertEqual(round(float(res["sinuosity"].max()), 3), 2.747)
