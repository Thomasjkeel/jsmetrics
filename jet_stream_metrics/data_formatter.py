# -*- coding: utf-8 -*-
"""
    Data formatter class for interacting, subsetting and calculating metrics from climate model outputs.

    Q: Why make a data formatter class and not just use functions to handle and format the data?
    A: When the data is in an object form, it can ... . Also, allows for information hiding.
"""

from . import compute_jetstream_metric

__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"
__description__ = "Data formatter class for interacting, subsetting and calculating metrics from climate model outputs"


class DataFormatter:
    """
    The DataFormatter object ...
    (see https://www.datacamp.com/community/tutorials/docstrings-python for docstring format)
    """
    def __init__(self, data):
        self.data = data
        self.variables = self.get_variable_list()

    @classmethod
    def with_available_metrics(self, data, all_metrics):
        self.data = data
        self.variables = self.get_variable_list()
        self.get_available_metrics(self, all_metrics)

    def get_available_metrics(self, all_metrics):
        self.all_metrics = all_metrics
        self.available_metrics = compute_jetstream_metric.get_available_metric_list(
            self.data, all_metrics)
        print("%s metrics available for this dataset:" %
              (len(self.available_metrics)))
        print("Metrics available:", self.available_metrics)

    def get_variable_list(self):
        variable_list = []
        for var in self.data.keys():
            if not '_bnds' in var:
                variable_list.append(var)
        return variable_list

    def subset(self, inplace=False, **kwargs):
        subset_data = self.data.copy()
        subset_data = subset_data.sel(**kwargs)
        if inplace:
            if hasattr(self, 'all_metrics'):
                return DataFormatter.with_available_metrics(subset_data, self.all_metrics)

        return DataFormatter(subset_data)


    def calc_metric(self):
        return
