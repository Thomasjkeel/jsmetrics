try:
    from jetstream_metrics import compute_jetstream_metric, data_formatter
except:
    from jetstream_metrics_recipe.jet_stream_metrics import compute_jetstream_metric, data_formatter
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

def main():
    print("Starting!")
    all_metrics = compute_jetstream_metric.JETSTREAM_METRICS
    data_dir = 'data/'
    UKESM1_SSP585_U = xr.open_dataset(data_dir + "ua_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    UKESM1_SSP585_V = xr.open_dataset(data_dir + "va_day_UKESM1-0-LL_ssp585_r2i1p1f2_gn_20150101-20491230.nc")
    UKESM1_SSP585 = xr.merge([UKESM1_SSP585_U, UKESM1_SSP585_V])
    ukesm1_ssp585 = data_formatter.DataFormatter(UKESM1_SSP585)
    ukesm1_ssp585 = ukesm1_ssp585.subset(lat=slice(0, 90), plev=25000)
    ukesm1_ssp585.get_available_metrics(all_metrics, return_coord_error=True)
    print("done!")
    return

