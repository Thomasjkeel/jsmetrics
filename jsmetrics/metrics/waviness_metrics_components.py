# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the waviness metrics used to
    characterise waviness in jet-stream in the literature.

    This file is in order of publish year of the metrics
    (see jsmetrics.details_for_all_metrics)

    Functions and Classes ordered by the year of the paper that uses each metric component first
"""

# imports
from jsmetrics.utils import spatial_utils


# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def calc_meridional_circulation_index(data):
    """
    Calculates the Meridional Circulation Index (MCI):
    MCI = v * abs(v) / u**2 * v**2
    When MCI = 0, the wind is purely zonal, and when MCI= 1 (-1), the flow is
    from the South (North).
    Component of method from Francis and Vavrus (2015) https://doi.org/10.1088/1748-9326/10/1/014005

    NOTE: The paper is not clear about whether the absolute value for MCI
    is taken instead thus 0-1

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u-component and v-component wind-speed values

    Returns
    ----------
    output : xarray.DataArray
        Array of Meridional Circulation Index (MCI) values
    """
    assert (
        "ua" in data.variables and "va" in data.variables
    ), "Cannot compute metric. 'ua' and/or 'va' not found in data"
    return data["va"] * abs(data["va"]) / (data["ua"] ** 2 + data["va"] ** 2)


def get_sinuosity_of_zonal_mean_zg(row, latitude_circle):
    """
    Works on a grouped data set and will calculate sinuosity of zonal mean
    geopotential (ZG) contour compared to a latitude circle

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    row : xarray.Dataset
        Data containing geopotential height (zg) and zonal_mean_zg_30Nto70N values
    latitude_circle : shapely.geometry.LineString
        Linestring of latitude circle around a hemisphere

    Returns
    ----------
    row : xarray.Dataset
        Data containing sinuousity value (determined by calc_great_circle_sinuosity function)
    """
    row["sinuosity"] = calc_great_circle_sinuosity(
        spatial_utils.get_one_contour_linestring(
            row["zg"], row["zonal_mean_zg_30Nto70N"].data
        ),
        latitude_circle,
    )
    return row


def calc_great_circle_sinuosity(line1, line2):
    """
    Calculates sinuosity by comparing the great circle distance between two (multi-)linestrings.

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    line1 : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Line to calculate great circle sinuosity from

    line2 : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Line to calculate great circle distance from

    Returns
    ----------
    sinuosity : int or float
        Sinuosity value between length of two lines
    """
    return spatial_utils.calc_total_great_circle_distance_along_line(
        line1
    ) / spatial_utils.calc_total_great_circle_distance_along_line(line2)
