# -*- coding: utf-8 -*-

"""
    General utility functions needed for the jet-stream metrics
    TODO: sort this out (maybe move to init file)
"""

### imports
import numpy as np
import itertools

### docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def make_climatology(data, freq):
    """
        Makes a climatology at given interval (i.e. days, months, season)
        
        Parameters
        ----------
        data (xarray.Dataset): data with regular time stamp
        freq (str): 'day', 'month' or 'season'
        
        Usage
        ----------
        climatology = make_climatology(data, 'month')
        
        
    """
    climatology = data.groupby("time.%s" % (freq)).mean("time")
    return climatology


def is_djf(month):
    """
        Mask used for getting DJF
    """
    return (month == 12) | (month >= 1) & (month <= 2)


def remove_duplicates(vals):
    """
        removes duplicates see: https://stackoverflow.com/questions/2213923/removing-duplicates-from-a-list-of-lists

        Used in a few metrics
    """

    vals.sort()
    vals = list(v for v,_ in itertools.groupby(vals))
    return vals


def get_all_plev_hPa(data):
    """
        Will get a list of all the pressure levels in the data in hPa 
    """
    if not 'plev' in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    plevs = np.array([plev for plev in data['plev']])
    if data['plev'].units == 'Pa':
        plevs = plevs/100 
    return plevs


def get_num_of_decimal_places(num):
    """
        func for getting number of decimal places in a float
        FOR GENERAL UTILS
    """
    num = '{:f}'.format(num).rstrip('0')
    decimal_places = num[::-1].find('.')
    if decimal_places < 0:
        decimal_places = 0
    return decimal_places