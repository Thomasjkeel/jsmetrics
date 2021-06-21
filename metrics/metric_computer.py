# -*- coding: utf-8 -*-
"""
    Metric Computer class for interacting, subsetting and calculating metrics from climate model outputs.

    Q: Why make a data formatter class and not just use functions to handle and format the data?
    A: When the data is in an object form, it can ... . Also, allows for information hiding.
"""

from . import compute_metrics

__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"
__description__ = "Data formatter class for interacting, subsetting and calculating metrics from climate model outputs"


class MetricComputer:
    """
    The MetricComputer object ...
    (see https://www.datacamp.com/community/tutorials/docstrings-python for docstring format)
    """
    def __init__(self, data):
        self.data = data
        self.get_variable_list()
        self.swap_all_coords()

    @classmethod
    def with_available_metrics(cls, data, all_metrics):
        data_with_available_metrics = cls(data)
        data_with_available_metrics.get_available_metrics(all_metrics)
        return data_with_available_metrics

    def get_available_metrics(self, all_metrics, return_coord_error=False):
        self.all_metrics = all_metrics
        self.available_metrics = compute_metrics.get_available_metric_list(
            self.data, all_metrics, return_coord_error)
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
            self.data = compute_metrics.swap_coord_order(self.data, coord)

    def subset(self, inplace=False, **kwargs):
        """
            Exposes the xarray .sel function
        """
        subset_data = self.data.copy()
        subset_data = subset_data.sel(**kwargs)
        if inplace:
            if hasattr(self, 'all_metrics'):
                return MetricComputer.with_available_metrics(subset_data, self.all_metrics)

        return MetricComputer(subset_data)

    def compute_metric_from_data(self, metric_name, **kwargs):
        """
            TODO: maybe add catch. Will print which metrics are available
        """
        result = compute_metrics.compute_metric(self.data, metric_name, **kwargs)
        return result
    
    def compute_all_metrics(self):
        """
            will go through and compute all metric which are available
        """
        if not hasattr(self, 'available_metrics'):
            print('please run .get_available_metrics() first')
            return
        