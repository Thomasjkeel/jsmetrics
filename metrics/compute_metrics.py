# -*- coding: utf-8 -*-

"""
    Contains the MetricComputer class and functions that run inside that class for subsetting
    and computing metrics from standardised netcdf data
"""

import xarray as xr
from .jetstream_metrics_dict import JETSTREAM_METRIC_DICT
from .general_utils import check_kwargs

__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


class MetricComputer:
    """
    Metric Computer class for interacting, subsetting and calculating metrics
    from climate model outputs.

    Q: Why make a data formatter class and not just use functions to handle and format the data?
    A: When the data is in an object form, it can ... . Also, allows for information hiding.

    (see https://www.datacamp.com/community/tutorials/docstrings-python for docstring format)
    """
    def __init__(self, data, all_metrics):
        assert isinstance(data, xr.Dataset), "data needs to be xarray.dataset"
        assert isinstance(all_metrics, dict) and len(all_metrics) > 0,\
                     "all metrics needs to be a dict with at least one value"
        self.data = data
        self.get_variable_list()
        self.swap_all_coords()
        self.all_metrics = all_metrics

    @classmethod
    def with_available_metrics(cls, data, all_metrics):
        data_with_available_metrics = cls(data, all_metrics)
        data_with_available_metrics.get_available_metrics(all_metrics)
        return data_with_available_metrics

    def get_available_metrics(self, return_coord_error=False):
        self.available_metrics = get_available_metric_list(
            self.data, self.all_metrics, return_coord_error)
        print("%s metrics available for this dataset:" %
              (len(self.available_metrics)))
        print("Metrics available:", self.available_metrics)

    def get_variable_list(self):
        self.variable_list = []
        for var in self.data.keys():
            if not '_bnds' in var:
                self.variable_list.append(var)

    def swap_all_coords(self):
        for coord in self.data.coords:
            if not self.data[coord].count() == 1:
                self.data = swap_coord_order(self.data, coord)

    def sel(self, **kwargs):
        """
            Exposes the xarray .sel function
        """
        new_data = self.data.copy()
        new_data = new_data.sel(**kwargs)
        if hasattr(self, 'available_metrics'):
            return MetricComputer.with_available_metrics(new_data, self.all_metrics)

        return MetricComputer(new_data, self.all_metrics)

    def isel(self, **kwargs):
        """
            Exposes the xarray .sel function
        """
        new_data = self.data.copy()
        new_data = new_data.isel(**kwargs)
        if hasattr(self, 'available_metrics'):
            return MetricComputer.with_available_metrics(new_data, self.all_metrics)
        else:
            return MetricComputer(new_data, self.all_metrics)

    def compute_metric_from_data(self, metric_name, calc_kwargs=None, subset_kwargs=None):
        """
            TODO: maybe add catch. Will print which metrics are available
        """
        calc_kwargs = check_kwargs(calc_kwargs)
        subset_kwargs = check_kwargs(subset_kwargs)
        if not hasattr(self, 'available_metrics'):
            try:
                self.get_available_metrics()
            except Exception as e:
                raise KeyError('A dictionary of all metrics is required') from e
        result = compute_metric(self.data, metric_name, all_metrics=self.all_metrics,\
                     subset_kwargs=subset_kwargs, calc_kwargs=calc_kwargs)
        return result

    def compute_all_metrics(self):
        """
            will go through and compute all metric which are available
        """
        if not hasattr(self, 'available_metrics'):
            try:
                self.get_available_metrics()
            except Exception as e:
                raise KeyError('A dictionary of all metrics is required.') from e


def subset_data(data, metric, ignore_coords=None):
    """
        Will subset the data based on the metric chosen
        TODO add way of only subsetting some coords data partially

        Parameters
        ----------
        data : xarray.Dataset
            climate data
        metric : dict
            jetstream metric from jetstream metric dictionary
        ignore_coords : list or tuple
            coordiantes to not subset
    """
    if ignore_coords:
        assert isinstance(ignore_coords, (list, tuple)), "ignore coords need to be list or set"
    assert isinstance(data, xr.Dataset), "data needs to be xarray.dataset"
    ## overwrite which coords will be changed
    coords_to_subset = get_coords_to_subset(ignore_coords, metric)
    ## check if subset is still possible
    if len(coords_to_subset) != 0:
        print('Subsetting data...')
        subset = data.copy()
        for coord in metric['coords'].keys():
            if coord in coords_to_subset:
                min_val = float(metric['coords'][coord][0])
                max_val = float(metric['coords'][coord][1])
                selection = {coord:slice(min_val, max_val)}
                subset = subset.sel(selection)
        subset = flatten_dims(subset)
        return subset
    else:
        return data


def get_coords_to_subset(ignore_coords, metric):
    if ignore_coords:
        coords_to_subset = set(metric['coords'].keys())
        for removed_coord in coords_to_subset.intersection(ignore_coords):
            print('Note:', removed_coord, 'has not been subset for the experiment')
        coords_to_subset = coords_to_subset.difference(set(ignore_coords))
        coords_to_subset = list(coords_to_subset)
    else:
        coords_to_subset = list(metric['coords'].keys())
    return coords_to_subset


def flatten_dims(data):
    """
        Supports subset and will flatten coordinates of an Xarray DataSet/DataArray
        with one value (so they are standardised)
    """
    for dim in data.dims:
        if data.dims[dim] == 1:
            selection = {dim:0}
            data = data.isel(selection)
    return data


def swap_coord_order(data, coord, ascending=True):
    """
        Will reverse the dimension if a higher number is first

        Parameters
        ----------
        data : xarray.Dataset
            climate data
        coord : str
            name from coord to change

        Useage
        ----------
        new_data = swap_coord_order(data, "lat")
    """
    first_val = 0
    last_val = -1
    if not ascending:
        first_val = -1
        last_val = 0
    if data[coord][first_val] > data[coord][last_val]:
        data = data.reindex(**{coord:list(reversed(data[coord]))})
    return data


def compute_metric(data, metric_name, all_metrics=None, return_coord_error=False,\
                 subset_kwargs=None, calc_kwargs=None):
    """
        Write function description

        Parameters
        ----------
        data : xarray.Dataset
            climate data
        metric_name : str
            name from jetstream metric file
    """
    assert isinstance(data, xr.Dataset), "data needs to be xarray.dataset"
    calc_kwargs = check_kwargs(calc_kwargs)
    subset_kwargs = check_kwargs(subset_kwargs)
    if not all_metrics:
        print('No metrics provided, defaulting to local JETSTREAM_METRICS file')
        all_metrics = JETSTREAM_METRIC_DICT
    assert isinstance(all_metrics, dict), "all metrics needs to be a dict with at least one value"
    ## check that you can actually compute metrics
    if check_all_coords_available(data, all_metrics[metric_name], return_coord_error)[0]\
                                 and check_all_variables_available(data, all_metrics[metric_name]):
        # print('all checks passed')
        pass
    else:
        print('cannot calculate %s metric from data provided' % (metric_name))
        # TODO have this return a useful message
        return

    ## subset data for metric
    subset = subset_data(data, all_metrics[metric_name], **subset_kwargs)

    ## calculate metric
    if subset:
        result = all_metrics[metric_name]['metric'](subset, **calc_kwargs)
    else:
        print('could not calculate metric from subsetted data')
        return False
    return result


def get_available_metric_list(data, all_metrics=None, return_coord_error=False):
    """
    Checks which variables can be used by the data

        Parameters
        ----------
        data : xr.Dataset or similar
            Xarray dataset

        all_metrics : dict (default: None)
            dictionary of jet-stream metrics

        return_coord_error : bool
            whether a message about where the correct coords but wrong
            coord values should be returned in available metrics list
            e.g. wrong pressure level (plev)

        Returns
        -------
        metric_available_list : list

        Usage
        -----
        m_list = get_available_metric_list(vwind_data, js_metrics)

    """
    assert isinstance(data, xr.Dataset), "data needs to be xarray.dataset"
    if not all_metrics:
        print('No metrics provided, defaulting to local JETSTREAM_METRICS file')
        all_metrics = JETSTREAM_METRIC_DICT

    available_metrics = []
    for metric_name in all_metrics:
        metric_is_usuable = {metric_name: 'usuable'}
        if check_all_variables_available(data, metric=all_metrics[metric_name]):
            # check that all coords exists in xarray data i.e. plev, lat, etc.
            metric_usable, coord_error_message = check_all_coords_available(data,\
                                                 all_metrics[metric_name], return_coord_error)
            ## will make return error message
            if return_coord_error and len(coord_error_message) > 0:
                metric_is_usuable =  {metric_name: "To use this metric" + coord_error_message}
            if metric_usable:
                available_metrics.append(metric_is_usuable)

    return available_metrics


def check_all_variables_available(data, metric):
    """
        Checks if all variables required to compute metric
        exist in the data.
    """
    for var in metric['variables']:
        if var in data.variables:
            pass
        else:
            return False
    return True


def check_all_coords_available(data, metric, return_coord_error=False):
    """
        Checks if all coords required to compute metric
        exist in the data.
    """
    assert isinstance(data, xr.Dataset), "data needs to be xarray.dataset"
    coord_error_message = ""
    metric_usable = True
    try: # TODO
        assert len(metric['coords']) >= 1, "Metric dictionary has less than 1 coordinate"
    except Exception as e:
        raise ValueError("Metric has no coordinates to subset. Usable: %s" % (metric_usable)) from e

    ## Loop over each coordinate in all metric dictionary and
    ## check if the coords exist in data and can be used for the metric calculation
    for coord in metric['coords'].keys():
        if coord in data.coords:
            coord_vals = metric['coords'][coord]
            coord_available = check_if_coord_vals_meet_reqs(data, coord, coord_vals)
            # if coord fails check, provide user information why
            if return_coord_error and not coord_available:
                coord_error_message += " '%s' needs to be between %s and %s." %\
                                     (str(coord), str(coord_vals[0]), str(coord_vals[1]))
            elif not coord_available:
                metric_usable = False
                break
        else:
            # if it does not exist then break loop as it is required for the metric
            metric_usable = False
            break
    return metric_usable, coord_error_message


def check_if_coord_vals_meet_reqs(data, coord, coord_vals):
    """
        Checks if the data has the correct coordinate values required
        for the metric.
    """
    min_val = float(coord_vals[0])
    max_val = float(coord_vals[1])
    ## check that the data is more than one value
    if data[coord].count() > 1:
        coord_val_avaialable = data[coord].loc[min_val: max_val]
        if len(coord_val_avaialable) == 0:
            return False

        return True
    else:
        return data[coord].values > min_val and data[coord].values < max_val
