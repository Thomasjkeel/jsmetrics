# -*- coding: utf-8 -*-

"""
    Plotting classes and functions
"""

# imports
import xarray as xr

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


class Plotter:
    """
    Base class for plotter. Input has to be xr.DataArray
    """

    def __init__(self, data):
        assert isinstance(
            data, xr.DataArray
        ), "input data to Plotter() class needs be xarray.DataArray type"
        self.data = data
        self._check_required_coords_are_in_data()

    def __init_subclass__(cls, req_coords, *a, **kw):
        cls.req_coords = req_coords

    def __repr__(self):
        return repr(self.data)

    def _check_required_coords_are_in_data(self):
        if not hasattr(self, "req_coords"):
            print("Using base Plotter() class")
            self.req_coords = ()
        for coord in self.req_coords:
            assert (
                coord in self.data.dims
            ), "'%s' not in data and is required for Plotter" % (coord)


class TimePlotter(Plotter, req_coords=("time",)):
    def time_grouper(self, time_group):
        """
        plotter.time_grouper('time.day').mean()
        """
        return self.data.groupby(time_group)

    def time_period_sel(self, starttime, endtime, inplace=False):
        try:
            if inplace:
                self.data = self.data.sel(time=slice(starttime, endtime))
                return self.data
            else:
                return self.data.sel(time=slice(starttime, endtime))
        except Exception:
            pass
            if inplace:
                self.data = self.data.isel(time=slice(starttime, endtime))
                return self.data
            else:
                return self.data.isel(time=slice(starttime, endtime))

    def time_period_mean(
        self, starttime, endtime, mean_dims=(...), inplace=False
    ):
        """
        time_plotter = TimePlotter(data_array)
        time_plotter.time_period_mean(starttime=0,\
            endtime=2, mean_dims=('time'), inplace=True)
        """
        if isinstance(starttime, int) or isinstance(starttime, float):
            starttime = self.data["time"].data[starttime]
            endtime = self.data["time"].data[endtime]

        if inplace:
            self.data = self.time_period_sel(starttime, endtime, inplace).mean(
                mean_dims
            )
            print("data updated")
            data = self.data
        else:
            data = self.time_period_sel(starttime, endtime).mean(mean_dims)
        data["time_selection"] = "%s to %s" % (starttime, endtime)
        return data

    def time_rolling_mean(self, time_periods, center=True):
        return self.data.rolling(time=time_periods, center=center).mean()

    def plot_time_series(self, **kwargs):
        """
        time_plotter.plot_time_series(ax=ax)
        """
        assert "time" in self.data.dims, "Time needs to be in data dimensions"
        assert self.data["time"].size > 1, "Data needs to have a least 2 dates"
        assert (
            len(self.data.dims) == 1
        ), "Data needs to have exactly 1 dimensions - time"
        return self.data.plot(**kwargs)

    def plot_rolling_mean(self, time_periods, center=True, **kwargs):
        """
        time_plotter.plot_time_series(ax=ax)
        time_plotter.plot_rolling_mean(ax=ax, time_periods=100)
        """
        rmean = self.time_rolling_mean(time_periods, center)
        return rmean.plot(**kwargs)


class TimeLatLonPlotter(TimePlotter, req_coords=("lat", "lon")):
    def map_data(self, **kwargs):
        """
        tll_plotter.time_period_mean(0,1000, mean_dims='time', inplace=True)
        tll_plotter.map()
        """
        if "time" in self.data.dims:
            assert (
                self.data["time"].size == 1
            ), "Data has too many dimensions.\
                Data needs to be lat and lon only"
        assert (
            "lat" in self.data.dims and "lon" in self.data.dims
        ), "Lat and Lon need to be in data"
        return self.data.plot(**kwargs)


class TimePlevLonPlotter(TimePlotter, req_coords=("plev", "lon")):
    pass
