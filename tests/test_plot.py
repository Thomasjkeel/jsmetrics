from metrics import compute_metrics
from metrics.jetstream_metrics_dict import JETSTREAM_METRIC_DICT
import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import unittest

class TestPlot(unittest.TestCase):
    def setUp(self):
        UKESM1_SSP585_U = xr.open_dataset("data/ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
        UKESM1_SSP585_V = xr.open_dataset("data/va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
        UKESM1_SSP585 = xr.merge([UKESM1_SSP585_U, UKESM1_SSP585_V])
        ukesm1_ssp585 = compute_metrics.MetricComputer(UKESM1_SSP585)
        self.data = ukesm1_ssp585.sel(lat=slice(0, 90))
        self.data = self.data.isel(time=slice(0,100))
        self.all_metrics = JETSTREAM_METRIC_DICT

    def plot_fig(self, save=False):
        max_lats = self.result[:,0]
        max_ws = self.result[:,1]
        fig, ax = plt.subplots(1)
        ax.plot(max_lats)
        ax.plot(max_ws)
        plt.legend(['latitude of max windspeed', 'windspeed'])
        if save:
            fig.savefig('experiments/figures/woolings_test.png', bbox_inches='tight')

    def test_plot(self): 
        self.result = self.data.compute_metric_from_data('Woolings2010', all_metrics=self.all_metrics, return_coord_error=False)
        self.assertIsInstance(self.result, np.ndarray)
        self.plot_fig()
        self.plot_fig(save=True)
        

if __name__ == "__main__":
    unittest.main()