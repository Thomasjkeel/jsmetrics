import xarray as xr
from metrics import compute_metrics
from metrics.jetstream_metrics_dict import JETSTREAM_METRIC_DICT


def load_uv_data(data_path):
    """
    Make this more general and add and test the data_dir
    """
    data_dir = "data/"
    UKESM1_SSP585_U = xr.open_dataset(
        data_dir + "ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc"
    )
    UKESM1_SSP585_V = xr.open_dataset(
        data_dir + "va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc"
    )
    UKESM1_SSP585 = xr.merge([UKESM1_SSP585_U, UKESM1_SSP585_V])
    ukesm1_ssp585 = compute_metrics.MetricComputer(
        UKESM1_SSP585, all_metrics=JETSTREAM_METRIC_DICT
    )
    ukesm1_ssp585 = ukesm1_ssp585.sel(lat=slice(0, 90))
    return ukesm1_ssp585


def main(data_path, **kwargs):
    print("Starting!")
    ukesm1_ssp585 = load_uv_data(data_path)
    ukesm1_ssp585.get_available_metrics(return_coord_error=True)
