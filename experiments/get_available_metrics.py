from numpy.lib.npyio import load
from metrics import compute_metrics, metric_computer
from metrics.jetstream_metrics_dict import JETSTREAM_METRIC_DICT
import xarray as xr


def load_uv_data(data_path):
    """
        Make this more general and add and test the data_dir
    """
    data_dir = 'data/'
    UKESM1_SSP585_U = xr.open_dataset(data_dir + "ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    UKESM1_SSP585_V = xr.open_dataset(data_dir + "va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    UKESM1_SSP585 = xr.merge([UKESM1_SSP585_U, UKESM1_SSP585_V])
    ukesm1_ssp585 = metric_computer.MetricComputer(UKESM1_SSP585)
    ukesm1_ssp585 = ukesm1_ssp585.subset(lat=slice(0, 90))
    return ukesm1_ssp585


def get_available_metrics(data, all_metrics):
    data.get_available_metrics(all_metrics, return_coord_error=True)
    

def main(data_path, **kwargs):
    print("Starting!")
    all_metrics = JETSTREAM_METRIC_DICT
    ukesm1_ssp585 = load_uv_data(data_path)
    get_available_metrics(ukesm1_ssp585, all_metrics)
    return 

