# -*- coding: utf-8 -*-

"""
    Tests for jet metrics (jetstream_metrics.py & jetstream_metrics_utils)
    Includes:
        – tests for each metric defined and expected outputs
        – test for each utility function used in the metric
"""

# imports
import unittest
from jsmetrics.metrics import jet_statistics, jet_statistics_components
from parameterized import parameterized
import xarray as xr
import numpy as np
from jsmetrics import (
    details_for_all_metrics,
)
from . import (
    set_up_test_uv_data,
    set_up_test_u_data,
    set_up_test_zg_data,
    set_up_nan_dataset,
)


# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


MAX_VARIABLES_IN_METRIC_DETAIL_DICT = 7


class TestMetricDetailsDict(unittest.TestCase):
    def setUp(self):
        self.metric_details = details_for_all_metrics.METRIC_DETAILS

    def test_metric_dict_keys(self):
        for metric_name in self.metric_details.keys():
            self.assertIsInstance(metric_name, str)

    def test_metric_dict_values(self):
        for metric in self.metric_details.values():
            self.assertIsInstance(metric, dict)
            self.assertEqual(len(metric.keys()), MAX_VARIABLES_IN_METRIC_DETAIL_DICT)
            self.assertListEqual(
                list(metric.keys()),
                [
                    "variables",
                    "coords",
                    "plev_units",
                    "metric",
                    "name",
                    "description",
                    "doi",
                ],
            )

    def test_variables(self):
        for metric in self.metric_details.values():
            self.assertIsInstance(metric["variables"], list)
            self.assertGreaterEqual(len(metric["variables"]), 0)
            self.assertLessEqual(
                len(metric["variables"]), MAX_VARIABLES_IN_METRIC_DETAIL_DICT
            )

    def test_metric_coords(self):
        for metric in self.metric_details.values():
            self.assertIsInstance(metric["coords"], dict)
            for coord in metric["coords"].keys():
                self.assertIsInstance(coord, str)

    @parameterized.expand(
        [
            ("plev", 0, 100000),
            ("lat", -91, 91),
            ("lon", -1, 361),
        ]
    )
    def test_each_coord(self, coord, min_value, max_value):
        for metric in self.metric_details.values():
            if coord in metric["coords"].keys():
                self.assertEqual(len(metric["coords"][coord]), 2)
                self.assertGreaterEqual(min(metric["coords"][coord]), min_value)
                self.assertLessEqual(max(metric["coords"][coord]), max_value)

    def test_funcs(self):
        for metric in self.metric_details.values():
            self.assertTrue(callable(metric["metric"]))


class TestMetricsOnOneDay(unittest.TestCase):
    def setUp(self):
        self.uvdata = set_up_test_uv_data()
        self.uvdata = self.uvdata.isel(time=0)
        self.zgdata = set_up_test_zg_data()
        self.zgdata = self.zgdata.isel(time=0)
        self.metric_details = details_for_all_metrics.METRIC_DETAILS

    def test_all_metrics(self):
        for metric_name in self.metric_details.keys():
            uvdata = self.uvdata.copy(deep=True)
            zgdata = self.zgdata.copy(deep=True)
            # do not include w10 or bp13 as they have a time window requirement
            if "Woollings2010" in metric_name or "BarnesPolvani2013" in metric_name:
                try:
                    self.metric_details[metric_name]["metric"](uvdata)
                except ValueError:
                    continue
            elif "Manney2011" in metric_name:
                self.metric_details[metric_name]["metric"](
                    uvdata, jet_core_plev_limit=(10000, 40000)
                )
            elif "Kuang2014" in metric_name or "Kerr2020" in metric_name:
                self.metric_details[metric_name]["metric"](uvdata.isel(plev=0))
                continue
            elif "zg" in self.metric_details[metric_name]["variables"]:
                self.metric_details[metric_name]["metric"](zgdata.isel(plev=0))
            else:
                self.metric_details[metric_name]["metric"](uvdata)


class TestArcherCaldeira2008(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        result = jet_statistics.archer_caldeira_2008(self.data)
        for col in [
            "mass_weighted_average_ws",
            "mass_flux_weighted_pressure",
            "mass_flux_weighted_latitude",
        ]:
            self.assertIn(col, result)
            # self.assertEqual(10, result[col].max())

        self.assertEqual(
            round(float(result["mass_weighted_average_ws"].max()), 4), 23.9048
        )

    def test_get_mass_weighted_average_windspeed(self):
        tested_func = jet_statistics.archer_caldeira_2008
        test_data = self.data.isel(plev=1)
        # should fail because needs two plevs
        self.assertRaises(ValueError, lambda: tested_func(test_data))
        self.assertRaises(
            KeyError,
            lambda: jet_statistics.archer_caldeira_2008(
                self.data.isel(time=0).drop_vars("time")
            ),
        )


class TestWoollings2010(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()
        # self.data = make_fake_seasonal_data(self.data)
        self.test_sig = self._get_test_sig()

    @staticmethod
    def _get_test_sig():
        # Seed the random number generator
        np.random.seed(42)
        time_step = 0.25
        period = 5.0
        time_vec = np.arange(0, 5, time_step)
        test_sig = np.sin(2 * np.pi / period * time_vec) + 0.5 * np.random.randn(
            time_vec.size
        )
        return test_sig

    def test_metric(self):
        tested_func = jet_statistics.woollings_et_al_2010
        result = tested_func(self.data, filter_freq=1, window_size=2)
        self.assertIsInstance(result, xr.Dataset)
        self.assertEqual(round(float(result["ff_jet_lat"][0]), 2), 36.25)
        self.assertEqual(round(float(result["ff_jet_speed"][0]), 2), 22.01)
        tested_func(self.data["ua"], filter_freq=1, window_size=2)

    def test_apply_lanczos_filter(self):
        tested_func = jet_statistics_components.apply_lanczos_filter
        test_ua = self.data["ua"]
        self.assertRaises(AssertionError, lambda: tested_func(self.data, 2, 4))
        self.assertRaises(ValueError, lambda: tested_func(test_ua, -2, 1))
        self.assertRaises(ValueError, lambda: tested_func(test_ua, 2, -1))
        self.assertRaises(ValueError, lambda: tested_func(test_ua, 2, 1))
        self.assertRaises(
            ValueError,
            lambda: tested_func(test_ua, test_ua["time"].count() + 2, 1),
        )
        self.assertRaises(
            ValueError,
            lambda: tested_func(test_ua, 2, test_ua["time"].count() + 1),
        )
        self.assertEqual(round(float(tested_func(test_ua, 2, 4).max()), 2), 99.18)

    def test_get_latitude_and_speed_where_max_ws(self):
        tested_func = jet_statistics_components.get_latitude_and_speed_where_max_ws
        self.assertRaises(AttributeError, lambda: tested_func(["lol"]))
        tested_data = self.data["ua"].isel(plev=0, lon=0, time=0)
        self.assertEqual(tested_func(tested_data)[0], 81.25)
        self.assertEqual(tested_func(tested_data)[1], -9.87109375)
        self.assertRaises(
            KeyError, lambda: tested_func(tested_data.rename({"lat": "lt"}))
        )
        nan_dataset = set_up_nan_dataset()
        self.assertEqual(tested_func(nan_dataset), (np.nan, np.nan))

    def test_apply_fourier_filter(self):
        tested_func = jet_statistics_components.apply_low_freq_fourier_filter
        res = tested_func(self.test_sig, 2)
        self.assertAlmostEqual(res[10].real, -0.0956296962962675, places=7)


class TestBarnesPolvani2013(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        result = jet_statistics.barnes_polvani_2013(
            self.data, filter_freq=1, window_size=2
        )
        self.assertIsInstance(result, xr.Dataset)
        self.assertEqual(round(float(result["jet_lat"][1]), 2), 35.97)
        self.assertEqual(round(float(result["jet_speed"][2]), 2), 16.82)
        self.assertEqual(float(result["jet_width"][1]), 17.5)

    def test_calc_jet_width_for_one_day(self):
        test_func = jet_statistics_components.calc_jet_width_for_one_day
        self.assertTrue(
            np.isnan(test_func(self.data["ua"].isel(time=0, plev=0), 25, None))
        )
        self.assertTrue(
            np.isnan(
                test_func(
                    self.data["ua"].isel(time=0, plev=0, lon=0, lat=slice(2, 5)),
                    0,
                    1,
                )
            )
        )

    def test_get_3_latitudes_and_speed_around_max_ws(self):
        test_func = jet_statistics_components.get_3_latitudes_and_speed_around_max_ws
        test_data = self.data["ua"].isel(time=0, plev=0, lat=slice(0, 2))
        res = test_func(test_data["lat"])
        self.assertEqual(len(res[0]), 3)
        self.assertTrue(np.isnan(res[0][-1]))
        test_data = self.data["ua"].isel(time=0, plev=0, lat=slice(68, 73))
        res = test_func(test_data["lat"])
        self.assertEqual(len(res[0]), 3)
        self.assertTrue(np.isnan(res[0][-1]))

    def test_get_3_neighbouring_coord_values(self):
        test_func = jet_statistics_components.get_3_neighbouring_coord_values
        res = test_func(45, 1)
        self.assertEqual(list(res), [44.0, 45.0, 46.0])


class TestScreenSimmonds2013(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_zg_data()

    def test_metric(self):
        # result = jetstream_metrics.screen_and_simmonds_2013(self.data)
        pass


class TestGrisePolvani2014(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        result = jet_statistics.grise_polvani_2014(self.data)
        jet_statistics.grise_polvani_2014(self.data["ua"])
        self.assertEqual(float(result["jet_lat"].min()), 35.38)
        self.assertEqual(float(result["jet_lat"].max()), 36.41)
        self.assertEqual(round(float(result["jet_speed"].max()), 5), 22.92644)
        self.assertRaises(
            KeyError,
            lambda: jet_statistics.grise_polvani_2014(
                self.data["ua"].isel(time=0).drop_vars("time")
            ),
        )


class TestBarnesPolvani2015(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        result = jet_statistics.barnes_polvani_2015(self.data)
        self.assertEqual(round(float(result["jet_lat"].mean()), 5), 43.19160)
        self.assertEqual(round(float(result["jet_speed"].max()), 5), 14.31844)
        self.assertEqual(round(float(result["jet_speed"].min()), 5), 13.52508)


class TestBarnesSimpson2017(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        test_func = jet_statistics.barnes_simpson_2017
        result = test_func(self.data)
        self.assertEqual(round(float(result["jet_lat"].mean()), 5), 36.25)
        self.assertEqual(round(float(result["jet_speed"].max()), 5), 22.05004)
        self.assertRaises(
            KeyError,
            lambda: jet_statistics.barnes_simpson_2017(
                self.data["ua"].isel(time=0).drop_vars("time")
            ),
        )


class TestBracegirdle2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        tested_func = jet_statistics.bracegirdle_et_al_2018
        test_data = self.data.sel(plev=slice(85000, 85000))
        result = tested_func(test_data)
        tested_func(test_data["ua"])
        # self.assertRaises(ValueError, lambda: tested_func(self.data))
        self.assertEqual(float(result["seasonal_JPOS"].max()), 37.725)
        self.assertEqual(float(result["annual_JPOS"].max()), 37.725)
        self.assertEqual(round(float(result["seasonal_JSTR"].max()), 3), 8.589)
        self.assertEqual(round(float(result["annual_JSTR"].max()), 3), 8.589)
        self.assertRaises(
            KeyError,
            lambda: jet_statistics.bracegirdle_et_al_2018(
                self.data["ua"].isel(time=0).drop_vars("time")
            ),
        )


class TestCeppi2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        tested_func = jet_statistics.ceppi_et_al_2018
        result = tested_func(self.data)
        self.assertEqual(float(result["jet_lat"][0].data), 37.316638365674194)
        self.assertEqual(float(result["jet_speed"][0].data), 22.341136932373047)

    def test_one_latlon_coord(self):
        tested_func = jet_statistics.ceppi_et_al_2018
        self.assertRaises(ValueError, lambda: tested_func(self.data.isel(lon=0)))
        self.assertRaises(ValueError, lambda: tested_func(self.data.isel(lat=0)))
        self.assertRaises(
            ValueError, lambda: tested_func(self.data.sel(lat=slice(0, 0)))
        )
        self.assertRaises(
            ValueError, lambda: tested_func(self.data.sel(lon=slice(0, 0)))
        )
        tested_func(self.data.sel(lon=slice(0, 0)), lon_resolution=1.875)
        tested_func(self.data.sel(lon=0), lon_resolution=1.875)
        self.assertRaises(
            ValueError,
            lambda: tested_func(
                self.data.sel(lon=slice(0, 0), lat=slice(0, 0)), lon_resolution=1.875
            ),
        )
        self.assertRaises(
            KeyError,
            lambda: jet_statistics.ceppi_et_al_2018(
                self.data["ua"].isel(time=0).drop_vars("time")
            ),
        )


class TestZappa2018(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_u_data()

    def test_metric(self):
        tested_func = jet_statistics.zappa_et_al_2018
        result = tested_func(self.data)
        self.assertEqual(round(float(result["jet_lat"][0].data), 3), 37.943)
        self.assertEqual(round(float(result["jet_speed"][0].data), 3), 22.341)

    def test_one_latlon_coord(self):
        tested_func = jet_statistics.zappa_et_al_2018
        self.assertRaises(ValueError, lambda: tested_func(self.data.isel(lon=0)))
        self.assertRaises(ValueError, lambda: tested_func(self.data.isel(lat=0)))
        self.assertRaises(
            ValueError, lambda: tested_func(self.data.sel(lat=slice(0, 0)))
        )
        self.assertRaises(
            ValueError, lambda: tested_func(self.data.sel(lon=slice(0, 0)))
        )
        tested_func(self.data.sel(lon=slice(0, 0)), lon_resolution=1.875)
        tested_func(self.data.sel(lon=0), lon_resolution=1.875)
        self.assertRaises(
            ValueError,
            lambda: tested_func(
                self.data.sel(lon=slice(0, 0), lat=slice(0, 0)), lon_resolution=1.875
            ),
        )
        self.assertRaises(
            KeyError,
            lambda: jet_statistics.zappa_et_al_2018(
                self.data["ua"].isel(time=0).drop_vars("time")
            ),
        )


class TestKerr2020(unittest.TestCase):
    def setUp(self):
        self.data = set_up_test_uv_data()

    def test_metric(self):
        tested_func = jet_statistics.kerr_et_al_2020
        test_data = self.data.sel(plev=50000)
        result = tested_func(test_data)
        self.assertEqual(result["jet_lat"].isel(time=0).dropna("lon").size, 192)
        self.assertEqual(float(result["jet_lat"].max()), 72.5)
        self.assertEqual(float(result["smoothed_jet_lat"].max()), 65.0)
        self.assertEqual(result["smoothed_jet_lat"].isel(time=0).dropna("lon").size, 61)
        self.assertRaises(
            KeyError,
            lambda: jet_statistics.kerr_et_al_2020(
                self.data["ua"].isel(time=0).drop_vars("time")
            ),
        )


if __name__ == "__main__":
    unittest.main()
