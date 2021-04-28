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
    # 3. 'coords' will contain the required CMIP6? standard coords and each coord
    #  will provide a list of 2 values: mininum value for coord and maximum
    #  value for coord and must be a number
        # 3.1 for 'plev' coord it is in millibars and higher pressure is
        #  considered the minimum value (e.g. 85000 mbar - 50000 mbar) 
    # 4. 'metric' is the name of a function in jetstream_metrics.py
JETSTREAM_METRICS = {"Woolings2010": {"variables": ["ua"], "coords": {"plev": [
    92500,  70000]}, "metric": jetstream_metrics.woolings_et_al_2010}}
# , "exact_coords": {"plev": [92500, 85000, 77500, 70000]}


def subset_data(data, metric):
    """
        Write function description
    """
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


def get_available_metric_list(data, all_metrics, return_coord_error=False):
    """
    Checks which variables can be used by the data
        
        Parameters
        ----------
        data : xr.Dataset or similar
            Xarray dataset 
        
        all_metrics : dict
            all jet-stream metrics

        return_coord_error : bool
            whether a message about where the correct coords but wrong
            coord values should be returned in available metrics list
            e.g. wrong pressure level (plev) 

        Returns
        -------
        metric_available_list : list
            
        
    """
    available_metrics = []
    for metric in all_metrics:
        metric_usable = True
        if check_all_variables_available(data, metric=all_metrics[metric]):
            # check that all coords exists in xarray data i.e. plev, lat, etc.
            if return_coord_error:
                coord_error_message = ""
            for coord in all_metrics[metric]['coords'].keys():
                if coord in data.coords:
                    coord_vals = all_metrics[metric]['coords'][coord]
                    coord_available = check_if_coord_vals_available(data, coord, coord_vals)
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
        else:
            metric_usable = False

        ## will make return error message
        if return_coord_error and len(coord_error_message) > 0:
            metric = metric + " – To use this metric" + coord_error_message
        if metric_usable:
            available_metrics.append(metric)
        

    return available_metrics


def get_usable_metrics_list(data, all_metrics):
    """
        Only looks at if correct variables exist for metrics and
        ignores lat/lon, temporal and spatial resolution and plev 
        TODO: add more informative 
    """
    usuable_metrics = []
    for metric in all_metrics:
        if check_all_variables_available(data, all_metrics[metric]):
            usuable_metrics.append(metric)
    return usuable_metrics


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


def check_if_coord_vals_available(data, coord, coord_vals):
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
        if coord == 'plev': ## TODO: think about changing this
            return data[coord].values > max_val and data[coord].values < min_val
        else:
            return data[coord].values > min_val and data[coord].values < max_val

    return True

