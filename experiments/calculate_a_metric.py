from metrics import compute_metrics, data_formatter
from metrics.jetstream_metrics_dict import JETSTREAM_METRIC_DICT
import xarray as xr


def main(data_path, metrics_to_use=None):
    print("Starting!")
    all_metrics = JETSTREAM_METRIC_DICT
    data_dir = 'data/'
    UKESM1_SSP585_U = xr.open_dataset(data_dir + "ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    UKESM1_SSP585_V = xr.open_dataset(data_dir + "va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    UKESM1_SSP585 = xr.merge([UKESM1_SSP585_U, UKESM1_SSP585_V])
    ukesm1_ssp585 = data_formatter.DataFormatter(UKESM1_SSP585)
    ukesm1_ssp585 = ukesm1_ssp585.subset(lat=slice(0, 90))
    # ukesm1_ssp585.get_available_metrics(all_metrics, return_coord_error=True)
    
    if metrics_to_use is None:
        print('Warning: No metric given. Using Woolings2010')
        metric_to_use = 'Woolings2010'
        result = ukesm1_ssp585.compute_metric_from_data(metric_to_use, all_metrics=all_metrics, return_coord_error=False)
    
    for  metric in metrics_to_use:
        print("calculating metric: %s" % (metric))
        try:
            # TODO: move this try and catch to the actual method
            result = ukesm1_ssp585.compute_metric_from_data(metric, all_metrics=all_metrics, return_coord_error=False)
            print('result:', result)
        except Exception as e:
            print("Unable to perform experiment. Error is:",e)
    print("done!")

