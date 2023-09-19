# -*- coding: utf-8 -*-
"""
    All spatial operations needed for the jet-stream metrics and algorithms

    Classes and Functions ordered alphabetically.
"""

import numpy as np
import xarray as xr
import collections
import math
import shapely.geometry
import matplotlib.pyplot


EARTH_RADIUS = 6371000.0  # m


def _guess_bounds(points, bound_position=0.5):
    """
    Guess bounds of grid cells.

    Simplified function from iris.coord.Coord.

    Author: Denis Sergev (https://github.com/dennissergeev) https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    Parameters
    ----------
    points: numpy.array
        Array of grid points of shape (N,).
    bound_position: float, optional
        Bounds offset relative to the grid cell centre.

    Returns
    -------
    Array of shape (N, 2).
    """
    diffs = np.diff(points)
    diffs = np.insert(diffs, 0, diffs[0])
    diffs = np.append(diffs, diffs[-1])

    diffs = _standardise_diffs_by_making_all_most_common_diff(diffs)

    min_bounds = points - diffs[:-1] * bound_position
    max_bounds = points + diffs[1:] * (1 - bound_position)

    return np.array([min_bounds, max_bounds]).transpose()


def _standardise_diffs_by_making_all_most_common_diff(diffs):
    """
    Lazy method to fill in gaps for bounds to make sure it is on a regular grid
    Adapted by githubuser:Thomasjkeel

    Author: Denis Sergev (https://github.com/dennissergeev) https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    """
    counter_of_diffs = collections.Counter(diffs)
    if len(counter_of_diffs) > 1:
        #  more than one difference found, so standardising to the most common
        print(
            "Warning: Standardising the guessed bounds of lat and lon to the most common bound width i.e. assumes regular grid (e.g. 1*1) and a gap in data"
        )
        most_common_key = list(collections.Counter(diffs).keys())[0]
        diffs = np.array([most_common_key] * len(diffs))
    return diffs


def _quadrant_area(radian_lat_bounds, radian_lon_bounds, radius_of_earth):
    """
    Calculate spherical segment areas.

    Taken from SciTools iris library.

    Area weights are calculated for each lat/lon cell as:
        .. math::
            r^2 (lon_1 - lon_0) ( sin(lat_1) - sin(lat_0))

    The resulting array will have a shape of
    *(radian_lat_bounds.shape[0], radian_lon_bounds.shape[0])*
    The calculations are done at 64 bit precision and the returned array
    will be of type numpy.float64.

    Author: Denis Sergev (https://github.com/dennissergeev) https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    Parameters
    ----------
    radian_lat_bounds: numpy.array
        Array of latitude bounds (radians) of shape (M, 2)
    radian_lon_bounds: numpy.array
        Array of longitude bounds (radians) of shape (N, 2)
    radius_of_earth: float
        Radius of the Earth (currently assumed spherical)

    Returns
    -------
    Array of grid cell areas of shape (M, N).
    """
    # ensure pairs of bounds
    if (
        radian_lat_bounds.shape[-1] != 2
        or radian_lon_bounds.shape[-1] != 2
        or radian_lat_bounds.ndim != 2
        or radian_lon_bounds.ndim != 2
    ):
        raise ValueError("Bounds must be [n,2] array")

    # fill in a new array of areas
    radius_sqr = radius_of_earth**2
    radian_lat_64 = radian_lat_bounds.astype(np.float64)
    radian_lon_64 = radian_lon_bounds.astype(np.float64)

    ylen = np.sin(radian_lat_64[:, 1]) - np.sin(radian_lat_64[:, 0])
    xlen = radian_lon_64[:, 1] - radian_lon_64[:, 0]
    areas = radius_sqr * np.outer(ylen, xlen)

    # we use abs because backwards bounds (min > max) give negative areas.
    return np.abs(areas)


def calc_spatial_integral(
    xr_da, lon_name="longitude", lat_name="latitude", radius=EARTH_RADIUS
):
    """
    Calculate spatial integral of xarray.DataArray with grid cell weighting.

    Author: Denis Sergev (https://github.com/dennissergeev) https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    Parameters
    ----------
    xr_da: xarray.DataArray
        Data to average
    lon_name: str, optional
        Name of x-coordinate
    lat_name: str, optional
        Name of y-coordinate
    radius: float
        Radius of the planet [metres], currently assumed spherical (not important anyway)

    Returns
    -------
    Spatially averaged xarray.DataArray.
    """
    lon = xr_da[lon_name].values
    lat = xr_da[lat_name].values

    area_weights = grid_cell_areas(lon, lat, radius=radius)

    return (xr_da * area_weights).sum(dim=[lon_name, lat_name])


def calc_spatial_mean(
    xr_da, lon_name="longitude", lat_name="latitude", radius=EARTH_RADIUS
):
    """
    Calculate spatial mean of xarray.DataArray with grid cell weighting.

    Author: Denis Sergev (https://github.com/dennissergeev) https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    Parameters
    ----------
    xr_da: xarray.DataArray
        Data to average
    lon_name: str, optional
        Name of x-coordinate
    lat_name: str, optional
        Name of y-coordinate
    radius: float
        Radius of the planet [metres], currently assumed spherical (not important anyway)

    Returns
    -------
    Spatially averaged xarray.DataArray.
    """
    lon = xr_da[lon_name].values
    lat = xr_da[lat_name].values

    area_weights = grid_cell_areas(lon, lat, radius=radius)
    aw_factor = area_weights / area_weights.max()

    return (xr_da * aw_factor).mean(dim=[lon_name, lat_name])


def calc_total_great_circle_distance_along_line(line):
    """
    Returns the total great circle (haversine) distance along a linestring
    or multilinestring

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    line : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Line to calculate great circle distance from

    Returns
    ----------
    total_distance : int or float
        Total distance in degrees or m

    """
    total_distance = 0
    if isinstance(line, shapely.geometry.multilinestring.MultiLineString):
        for i, _ in enumerate(line.geoms):
            total_distance += get_great_circle_distance_along_linestring(
                shapely.geometry.LineString((line.geoms[i]))
            )
    elif isinstance(line, shapely.geometry.LineString):
        total_distance += get_great_circle_distance_along_linestring(line)
    else:
        return np.nan
    return total_distance


def get_great_circle_distance_along_linestring(line):
    """
    Calculate great circle distance along the length of linestring

    Parameters
    ----------
    line : shapely.geometry.LineString
        Line to calculate great circle (haversine) distance along

    Returns
    ----------
    distance : float
        Great circle (haversine) distance along input line
    """
    numCoords = len(line.coords) - 1
    distance = 0
    for i in range(0, numCoords):
        point1 = line.coords[i]
        point2 = line.coords[i + 1]
        distance += haversine(point1[0], point1[1], point2[0], point2[1])
    return distance


def get_latitude_circle_linestring(latitude, lon_min, lon_max):
    """
    Will return a linestring of a latitude circle.

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    latitude : int or float
        given latitude to calculate circle from
    lon_min : int or float
        Minimum longitude for circle to extend to
    lon_max : int or float
        Maximum longitude for circle to extend to

    Returns
    ----------
    circle : shapely.geometry.LineString
        Linestring of latitude circle around a hemisphere
    """
    vals = np.column_stack(
        (
            np.arange(lon_min, lon_max + 0.1, 0.5),
            np.array([latitude] * len(np.arange(lon_min, lon_max + 0.1, 0.5))),
        )
    )
    circle = shapely.geometry.LineString(vals)
    return circle


def get_one_contour_linestring(dataarray, contour_level):
    """
    Returns a linestring or multi-linestring of a given contour.

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    dataarray : xarray.DataArray
        Array of Geopotential height (zg) values to calculate contour from
    contour_level :
        Value with which to calculate a contour from geopotential height

    Returns
    ----------
    contour_line : shapely.geometry.LineString or shapely.geometry.MultiLineString
        Contour line of geopotential height (zg) a given contour
    """
    assert isinstance(dataarray, xr.DataArray), "Data needs to be type xr.DataArray"
    assert (
        "lat" in dataarray.coords and "lon" in dataarray.coords
    ), "Data array needs to have latitude and longitude coords"
    one_contour = dataarray.plot.contour(levels=[contour_level])
    matplotlib.pyplot.close()
    one_contour_segments = seperate_one_contour_into_line_segments(
        one_contour.get_paths()[0]
    )
    contour_line = shapely.geometry.MultiLineString(one_contour_segments)
    return contour_line


def grid_cell_areas(lon1d, lat1d, radius=EARTH_RADIUS):
    """
    Calculate grid cell areas given 1D arrays of longitudes and latitudes
    for a planet with the given radius.

    Author: Denis Sergev (https://github.com/dennissergeev) https://gist.github.com/dennissergeev/60bf7b03443f1b2c8eb96ce0b1880150

    Parameters
    ----------
    lon1d: numpy.array
        Array of longitude points [degrees] of shape (M,)
    lat1d: numpy.array
        Array of latitude points [degrees] of shape (M,)
    radius: float, optional
        Radius of the planet [metres] (currently assumed spherical)

    Returns
    -------
    Array of grid cell areas [metres**2] of shape (M, N).
    """
    lon_bounds_radian = np.deg2rad(_guess_bounds(lon1d))
    lat_bounds_radian = np.deg2rad(_guess_bounds(lat1d))
    area = _quadrant_area(lat_bounds_radian, lon_bounds_radian, radius)
    return area


def haversine(lon1, lat1, lon2, lat2):
    """
    Calculate the great circle distance between two points
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (
        math.sin(dlat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(dlon / 2) ** 2
    )
    c = 2 * math.asin(math.sqrt(a))
    r = 6371  # Radius of earth in kilometers. Use 3956 for miles
    return c * r


def seperate_one_contour_into_line_segments(one_contour):
    """
    Seperates a list of vertices of a given contour into multiple segments for
    eventual conversion to multi-linestring.

    Needed with new Matplotlib depreciation to collections and allsegs of the contour plots.

    Component of method from Cattiaux et al (2016) https://doi.org/10.1002/2016GL070309

    Parameters
    ----------
    one_contour : matplotlib.path.Path
        Path object of one contour from a contour plot

    Returns
    ----------
    all_segments : list of lists
        List of all segments that can be converted to a multilinestring
    """

    all_segments = []
    current_segment = []
    for ind, (segment, code) in enumerate(one_contour.iter_segments()):
        if code == 1:
            if ind != 0:
                all_segments.append(current_segment)
            current_segment = []
        current_segment.append(segment)
    all_segments.append(current_segment)
    return all_segments
