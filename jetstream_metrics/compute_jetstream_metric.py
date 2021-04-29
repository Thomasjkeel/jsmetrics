# -*- coding: utf-8 -*-

"""
    Functions for subsetting and compute metrics from standardised climate model output data
"""

from . import jetstream_metrics

__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


# RULES for JETSTREAM_METRICS
    # 1. must have the keys: 'variables', 'coords' and 'metric'
    # 2. 'variables' will contain the required CMIP6? model output variable names
    # 3. 'coords' will contain the required CMIP6? standard coords and each coord will provide a list of 2 values:
    #  mininum value for coord and maximum value for coord and must be a number
        # 3.1 for 'plev' coord it is in millibars and higher pressure is
        #  considered the minimum value (e.g. 85000 mbar - 50000 mbar) 
    # 4. 'metric' is the name of a function in jetstream_metrics.py
    
JETSTREAM_METRICS = {"Woolings2010": {"variables": ["ua"], "coords": {"plev": [
    92500,  70000]}, "metric": jetstream_metrics.woolings_et_al_2010,
    "description":"Woolings et al. 2010 TODO"}} # , "exact_coords": {"plev": [92500, 85000, 77500, 70000]}

def subset_data(data, metric):
    """
        Will subset the data based on the metric chosen

        Parameters
        ----------
        data : xarray.Dataset
    """
    ## 
    # for coord in all_metrics[metric]['coords'].keys():
    #     pass
    # i.e. data.sel(plev=())
    # i.e. if the term is not found, then user input to slightly adjust or skip?
    
    # for key in JETSTREAM_METRICS[metric].keys():
    #  subset
    return


def compute_metric(data, metric):
    """
        Write function description
    """
    # subset_data(data, metric)
    # then
    # if JETSTREAM_METRICS[metric]["metric"]:
    # JETSTREAM_METRICS[metric]["metric"](data)
    return


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
    if not all_metrics:
        print('No metrics provided, defaulting to local JETSTREAM_METRICS file')
        all_metrics = JETSTREAM_METRICS

    available_metrics = []
    for metric_name in all_metrics:
        if check_all_variables_available(data, metric=all_metrics[metric_name]):
            # check that all coords exists in xarray data i.e. plev, lat, etc.
            metric_usable, coord_error_message = check_all_coords_available(data, all_metrics[metric_name], return_coord_error)

        ## will make return error message
        if return_coord_error and len(coord_error_message) > 0:
            metric_name = metric_name + " – To use this metric" + coord_error_message
        if metric_usable:
            available_metrics.append(metric_name)

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
    coord_error_message = ""
    metric_usable = True
    assert len(metric['coords']) >= 1, "Metric dictionary has less than 1 coordinate" # TODO

    ## Loop over each coordinate in all metric dictionary and check if the coords exist in data and can be used for the metric calculation
    for coord in metric['coords'].keys():
        if coord in data.coords:
            coord_vals = metric['coords'][coord]
            coord_available = check_if_coord_vals_meet_reqs(data, coord, coord_vals)
            # if coord fails check, provide user information why
            if return_coord_error and not coord_available:
                    coord_error_message += " the coord: %s needs to be between %s and %s." % (str(coord), str(coord_vals[0]), str(coord_vals[1]))
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
    else:
        if coord == 'plev' and min_val > max_val:
            return data[coord].values > max_val and data[coord].values < min_val
        else:
            return data[coord].values > min_val and data[coord].values < max_val

    return True
