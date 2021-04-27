# -*- coding: utf-8 -*-

"""
    Functions for subsetting and compute metrics from standardised climate model output data
"""

from . import jetstream_metrics

__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


# JETSTREAM_METRICS will have list of all metrics and how the data is being subset
JETSTREAM_METRICS = {"Woolings2010": {"variable": ["ua"], "plev": [
    "85000",  "70000"], "metric": jetstream_metrics.woolings_et_al_2010}}


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


def get_available_metric_list(data, all_metrics):
    available_metrics = []
    for metric in all_metrics:
        metric_usable = False
        for metric_property in all_metrics[metric].keys():
            # check that data exists in the data loaded by the class
            # if not metric_property in data.coord:
            # metric_usuable = False
            # break # TODO: check logic and test
            if metric_property in data.coords:
                print('this',metric_property)
            pass
        if metric_usable:
            available_metrics.append(metric)
    return available_metrics


