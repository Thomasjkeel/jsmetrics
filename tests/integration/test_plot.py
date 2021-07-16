# -*- coding: utf-8 -*-

"""
    Tests for plotting outputs from metrics
"""

### imports
from metrics import compute_metrics
from metrics.jetstream_metrics_dict import JETSTREAM_METRIC_DICT
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import unittest

class TestPlot(unittest.TestCase):
    def setUp(self):
        UKESM1_SSP585_U = xr.open_dataset("tests/data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20150105.nc")
        UKESM1_SSP585_V = xr.open_dataset("tests/data/va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20150105.nc")
        UKESM1_SSP585 = xr.merge([UKESM1_SSP585_U, UKESM1_SSP585_V])
        ## make fake seasonal data
        UKESM1_SSP585['time'] = np.array(['2015-01-01T00:00:00.000000000', '2015-01-02T00:00:00.000000000',
       '2015-02-01T00:00:00.000000000', '2015-02-02T00:00:00.000000000',
       '2015-02-03T00:00:00.000000000'], dtype='datetime64[ns]')
        self.metric_computer = compute_metrics.MetricComputer(UKESM1_SSP585, all_metrics=JETSTREAM_METRIC_DICT)

    def plot_fig(self, save=False):
        max_lats = self.result['max_lats']
        max_ws = self.result['max_ws']
        ff_max_lats = self.result['ff_max_lats']
        ff_max_ws = self.result['ff_max_ws']
        fig, ax = plt.subplots(1)
        ax.plot(max_lats)
        ax.plot(ff_max_lats)
        ax.plot(max_ws)
        ax.plot(ff_max_ws)
        plt.legend(['latitude', 'fourier filtered latitude', 'windspeed', 'fourier filtered wind-speed'])
        if save:
            fig.savefig('tests/figures/woolings_test.png', bbox_inches='tight')

    def test_plot(self): 
        self.result = self.metric_computer.compute_metric_from_data('Woolings2010', calc_kwargs={"filter_freq":1, "window_size":2})
        self.assertIsInstance(self.result, xr.Dataset)
        self.plot_fig()
        self.plot_fig(save=False)
        

if __name__ == "__main__":
    unittest.main()