# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to
    identify or classify jet-stream in the literature.

    This file is in order of publish year of the metrics
    (see details_for_all_metrics.py)

    Functions and Classes ordered by the year of the paper that uses each algorithm first.
"""

# imports
from jsmetrics.utils import data_utils, windspeed_utils
import numpy as np
import scipy.ndimage
import xarray as xr

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_sum_weighted_ws(data, all_plevs_hPa):
    r"""
    Get sum of weighted windspeed.
    sum weighted windspeed is calculated as follows:

    .. math::
        \int_{p1}^{p2} (u^2+v^2)^{1/2} \,dp

    where p1, p2 is min, max pressure level

    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev' and 'time'.
    all_plevs_hPa : array-like
        list of hPa unit pressure levels

    Returns
    ----------
    sum_weighted_ws : xarray.Dataset
        Data containing sum weighted windspeed values
    """
    if "plev" not in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError("array of pressure level needs to be list or numpy.array")

    sum_weighted_ws = 0
    for plev, (i, plev_hPa) in zip(data["plev"], enumerate(all_plevs_hPa)):
        if i != 0:
            plev_hPa = plev_hPa - all_plevs_hPa[i - 1]
        sum_weighted_ws += (
            (
                data.sel(plev=plev, method="nearest")["ua"] ** 2
                + data.sel(plev=plev, method="nearest")["va"] ** 2
            )
            ** (1 / 2)
        ) * plev_hPa
    return sum_weighted_ws


def get_weighted_average_ws(sum_weighted_ws, all_plevs_hPa):
    r"""
    Calculates weighted average windspeed as follows:

    .. math::
        \alpha vel =  \frac{1}{(p2-p1)} \int_{p1}^{p2} (u^2+v^2)^{1/2} \,dp

    where p1, p2 is min, max pressure level.

    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    Parameters
    ----------
    sum_weighted_ws : xarray.Dataset
        Data containing sum weighted windspeed values
    all_plevs_hPa : array-like
        list of hPa unit pressure levels

    Returns
    ----------
    weighted_average_ws : xarray.Dataset
        Data containing weighted average windspeed values
    """
    if not isinstance(all_plevs_hPa, (np.ndarray)):
        raise TypeError("array of pressure level needs to be a list or numpy.array")
    weighted_average_ws = sum_weighted_ws * (
        1 / (all_plevs_hPa.max() - all_plevs_hPa.min())
    )
    return weighted_average_ws


def get_all_hPa_list(data):
    """
    Will get a list of all the pressure levels in the data in hPa/mbar

    Parameters
    ----------
    data : xarray.Dataset
        data with plev coord

    Returns
    ----------
    plev : np.array
        arr containing pressure level list in data in hPa/mbar
    """
    if "plev" not in data.coords:
        raise KeyError("Data does not contain coord: 'plev'")

    if (
        data["plev"].units != "Pa"
        and data["plev"].units != "hPa"
        and data["plev"].units != "mbar"
        and data["plev"].units != "millibars"
    ):
        raise ValueError("Plev units need to be mbar, millibars, Pa or hPa")

    plevs = np.array([plev for plev in data["plev"]], dtype="float")
    if data["plev"].units == "Pa":
        plevs = plevs / 100
    return plevs


def get_local_jet_occurence(row, ws_threshold, u_threshold):
    r"""
    Each jet occurence is detected based on three rules applied to inputted
    wind speed (V = [u, v]):
        1. \|V\| is a local maxima in latitude and altitude plane
        2. \|V\| ≥ 30 m/s
        3. \|u\| ≥ 0 m/s.

    Component of method from Schiemann et al 2009 https://doi.org/10.1175/2008JCLI2625.1

    NOTE: will only work if 1 day is the resolution

    Parameters
    ----------
    row : xarray.Dataset
        Data of a single time unit containing windspeed (ws), plev, lat, lon
    ws_threshold : int or float
        Windspeed threshold used to extract jet events (default: 30 ms-1)
    u_threshold : int or float
        Windspeed threshold used to extract u-component wind speed (default: 0 ms^{-1})

    Returns
    ----------
    row : xarray.Dataset
        Data of a single time unit with value for jet-maxima (1 == maxima, 0 == none)

    """
    # Step 0. Squeeze time for method
    if "time" in row.dims:
        row = row.squeeze("time")

    row["jet_occurence"] = (
        ("plev", "lat", "lon"),
        np.zeros((row["plev"].size, row["lat"].size, row["lon"].size)),
    )
    all_jet_occurences = []
    for lon in row["lon"]:
        current = row.sel(lon=lon, method="nearest")

        # Finding local maxima indices in the 2D array
        maxima_indices_ax0 = np.column_stack(
            data_utils.get_local_maxima(current["ws"].data, axis=0)
        )
        maxima_indices_ax1 = np.column_stack(
            data_utils.get_local_maxima(current["ws"].data, axis=1)
        )
        # find intersection of the two arrays intersection
        maxima_indices = data_utils.find_intersection_between_two_array_of_arrays(
            maxima_indices_ax0, maxima_indices_ax1
        )

        # Filter indices to remove maximas that neighbour each other (taking the first instance)
        filtered_maxima_indices = data_utils.filter_local_extremes_to_min_distance(
            maxima_indices, min_distance_threshold=2
        )
        if len(filtered_maxima_indices) > 0:
            # Creating a boolean mask for the filtered maxima
            maxima_mask = np.zeros_like(current["ws"], dtype=bool)
            maxima_mask[
                filtered_maxima_indices[:, 0], filtered_maxima_indices[:, 1]
            ] = True

            # Mask coordinates of filtered maxima
            current["jet_occurence"] = current["ws"].where(
                (maxima_mask)
                & (abs(current["ws"]) >= ws_threshold)
                & (current["ua"] > u_threshold)
            )
        else:
            # Set all values to np.nan
            current["jet_occurence"] = current["jet_occurence"].where(
                current["jet_occurence"] > 0
            )
        all_jet_occurences.append(current["jet_occurence"])
    all_jet_occurences_dataarray = xr.concat(all_jet_occurences, dim="lon")

    # Set all non-NaN values to 1
    all_jet_occurences_dataarray = (
        all_jet_occurences_dataarray / all_jet_occurences_dataarray
    )
    row["jet_occurence"] = all_jet_occurences_dataarray.transpose("plev", "lat", "lon")

    return row


def run_jet_core_and_region_algorithm_on_one_day(
    row,
    jet_core_plev_limit,
    jet_core_ws_threshold,
    jet_boundary_ws_threshold,
    ws_drop_threshold,
    jet_core_lat_distance,
    check_diagonals,
):
    r"""
    This method detects jet cores and defines a boundary region beside those cores based on two windspeed thresholds.
    Two additional checks are applied after initial detection of cores to check whether cores within the same windspeed region
    are part of the same feature (see 'jet_boundary_ws_threshold'). This is achieved by checking whether regions with multiple
    jet cores are more than a certain distance apart (see 'jet_core_lat_distance') and the windspeed between two cores does
    not drop below a threshold (see 'ws_drop_threshold'). This function runs this method on a single time unit.

    This method returns four outputs
        1. 'jet_core_mask' -- Regions within each latitude/altitude that are local maxima have windspeeds above the 'jet_core_ws_threshold'
        2. 'jet_region_mask' -- Regions above, below, left and right of the jet core with windspeed above the 'jet_boundary_ws_threshold'
        3. 'jet_region_above_ws_threshold_mask' -- All contigious regions of windspeeds emcompassing a jet core above the 'jet_boundary_ws_threshold' (i.e. not just above, below, left and right)
        4. 'ws' -- Wind speed calculated from 'ua', 'va' inputs.

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    row : xarray.Dataset
        Data of single time unit containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev'.
    jet_core_plev_limit: tuple or array
        Sequence of two values relating to the pressure level limit of the jet cores (original paper uses (100, 400) hPa).
    jet_core_ws_threshold : int or float
        Threshold used for jet-stream core point.
    jet_boundary_ws_threshold : int or float
        Threshold for jet-stream boundary point.
    ws_drop_threshold : int or float
        Threshold for drop in windspeed along the line between cores.
    jet_core_lat_distance : int or float
        Threshold for maximum distance between cores to be counted the same.
    check_diagonals : bool
        Whether to check the diagonal edges of each maxima in the latitude-altitude plane.

    Returns
    ----------
    row : xarray.Dataset
        Data for one time unit containing four new variables (ws, jet_core_mask, jet_region_mask, jet_region_above_ws_threshold_mask)
    """
    # Step 0. Squeeze time for method
    if "time" in row.dims:
        row = row.squeeze("time")

    # Step 1. Get potential cores (to later subset) using the wind speed thresholds and jet core pressure level limit.
    row["potential_jet_cores"] = row["ws"].where(
        lambda val: (val >= jet_core_ws_threshold)
        & (val.plev >= min(jet_core_plev_limit))
        & (val.plev <= max(jet_core_plev_limit)),
        0,
    )
    row["potential_jet_regions"] = row["ws"].where(
        lambda val: val >= jet_boundary_ws_threshold, 0
    )

    # Step 2. Get local maximas at each longitude
    local_maximas_dict = get_local_maximas_at_each_longitude(
        row, check_diagonals, var_name="potential_jet_cores"
    )

    # Step 3. Loop through the local maximas and make an mask of initial jet cores (not official as some may be in same region)
    mask_shape = row["ws"].isel(lon=0).shape
    all_jet_core_mask = get_all_jet_core_mask(
        local_maximas_dict=local_maximas_dict, mask_shape=mask_shape
    )

    # Step 4. Loop through contigious regions (here known as jet region contours) and check if a jet core is within them
    all_jet_region_contour_mask = get_all_jet_region_contour_mask(
        row, local_maximas_dict=local_maximas_dict
    )

    # Step 5. Get the jet region around each local maxima only including above, below, left and right of maxima.
    all_jet_regions_mask = get_all_jet_regions_mask(
        row,
        all_jet_region_contour_mask=all_jet_region_contour_mask,
        local_maximas_dict=local_maximas_dict,
    )

    # Step 6. Remove old jet regions and define two outputs (jet_region, region_above_ws_threshold)
    row = row.drop_vars("potential_jet_regions")
    row["jet_region_mask"] = (
        ("plev", "lat", "lon"),
        np.clip(all_jet_regions_mask, 0, 1),
    )
    row["jet_region_contour_mask"] = (
        ("plev", "lat", "lon"),
        all_jet_region_contour_mask,
    )

    # Step 7. Run checks on jet cores, to check if they are part of the same jet feature
    jet_core_masks = run_checks_on_jet_cores_and_return_jet_cores(
        row,
        all_jet_core_mask,
        local_maximas_dict,
        jet_core_lat_distance,
        ws_drop_threshold,
    )

    # Step 8. Remove old and add actual jet core mask
    row = row.drop_vars("potential_jet_cores")
    row["jet_core_mask"] = (("plev", "lat", "lon"), jet_core_masks)
    return row


def get_local_maximas_at_each_longitude(row, check_diagonals, var_name):
    """
    Runs 'find_local_maxima_in_2d_dataarray' on each longitude (see docs for find_local_maxima_in_2d_dataarray for more information)

    Parameters
    ----------
    row : xr.DataArray
        Data of single time unit containing the variables: `var_name`, and the coordinates: 'lon', 'lat', 'plev'
    check_diagonals : bool
        Whether to check the diagonal edges of each maxima in the latitude-altitude plane.
    var_name : str
        Name of data array variable to find local maxima from e.g. 'potential_jet_cores'

    Returns
    ----------
    local_maximas_dict : dict
        Keys are longitude coordinates, values are local maxima locations

    """
    local_maximas_dict = {}
    for _, lon in enumerate(row["lon"]):
        pot_core_one_lon = row[var_name].sel(lon=lon)
        if check_diagonals:
            local_maximas = find_local_maxima_in_2d_dataarray_with_diagonals(
                pot_core_one_lon
            )
        else:
            local_maximas = find_local_maxima_in_2d_dataarray(pot_core_one_lon)
        local_maximas_dict[float(lon)] = local_maximas
    return local_maximas_dict


def find_local_maxima_in_2d_dataarray(arr):
    """
    Find indices of local maximas within a 2-D array. Should return two values which relate to the position of local maximas (if any).

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    arr : xr.DataArray
        A 2-D array with numeric values (i.e. dtypes float or int)

    Returns
    ----------
    local_maxima : xr.DataArray
        A 2-D array with values relating to the index of the local maximas

    Examples
    ----------
    .. code-block:: python

        import numpy as np
        import xarray as xr

        # Example array
        data = np.array([
            [1, 2, 3, 4, 5],
            [2, 3, 4, 5, 4],
            [3, 4, 5, 6, 7],
            [4, 5, 6, 7, 6],
            [5, 6, 7, 8, 9]
        ])

        # Convert NumPy array to xarray DataArray
        data_array = xr.DataArray(data)

        local_maxima_indices = find_local_maxima_in_2d_dataarray(data_array)

        for i, j in local_maxima_indices:
            print(f"Local maximum at ({i}, {j}): {data_array[i, j]}")

    """
    # Add padding to array to allow edges to be picked up
    arr = data_utils.add_pad_to_array(arr)

    # Calculate neighbors for all interior points
    neighbors = np.stack(
        [
            arr[:-2, 1:-1],  # Above
            arr[2:, 1:-1],  # Below
            arr[1:-1, :-2],  # Left
            arr[1:-1, 2:],  # Right
        ],
        axis=-1,
    )

    # Find local maximas
    interior_maxima = arr[1:-1, 1:-1] > np.max(neighbors, axis=-1)
    interior_indices = np.transpose(np.where(interior_maxima))

    return np.array(interior_indices)


def find_local_maxima_in_2d_dataarray_with_diagonals(arr, threshold=10):
    """
    Find indices of local maximas within a 2-D array with check for the diagonals.
    Should return two values which relate to the position of local maximas (if any).

    Component of method from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
    & Kuang et al (2014) (https://doi.org/10.1007/s00704-013-0994-x)

    Parameters
    ----------
    arr : xr.DataArray
        A 2-D array with numeric values (i.e. dtypes float or int)

    Returns
    ----------
    local_maxima : xr.DataArray
        A 1-D array with values relating to the index of the local maximas

    Examples
    ----------
    .. code-block:: python

        import numpy as np
        import xarray as xr

        # Example array
        data = np.array([
            [1, 2, 3, 4, 5],
            [2, 3, 4, 5, 4],
            [3, 4, 5, 6, 7],
            [4, 5, 6, 7, 6],
            [5, 6, 7, 8, 9]
        ])

        # Convert NumPy array to xarray DataArray
        data_array = xr.DataArray(data)

        local_maxima_indices = find_local_maxima_in_2d_dataarray_with_diagonals(data_array)

        for i, j in local_maxima_indices:
            print(f"Local maximum at ({i}, {j}): {data_array[i, j]}")

    """
    local_maxima = []

    # Add padding to array to allow edges to be picked up
    arr = data_utils.add_pad_to_array(arr)

    # Find local maximum values in a window around each element
    neighborhood_size = 3
    local_max = scipy.ndimage.maximum_filter(
        arr, footprint=np.ones((neighborhood_size, neighborhood_size))
    )

    # Compare local maximum with the threshold and original array
    local_maxima = (arr == local_max) & (arr > threshold)

    # Remove duplicates (this strategy can produce duplicates)
    local_maxima = np.unique(np.array(local_maxima), axis=0)

    # get index of local maxima
    local_maxima = np.dstack(np.where(local_maxima))[0] - 1
    return local_maxima


def get_all_jet_core_mask(local_maximas_dict, mask_shape):
    """
    Runs 'create_mask_using_local_maximas' using each

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    local_maximas_dict : dict
        Keys are longitude coordinates, values are local maxima locations

    mask_shape : tuple
        Values to create a mask of given shape from
        .
    Returns
    ----------
    all_jet_core_mask : numpy.array
        Mask of the jet cores

    """
    all_jet_core_mask = np.array([])
    for ind, local_maximas in enumerate(local_maximas_dict.values()):
        current_jet_core_mask = create_mask_using_local_maximas(
            local_maximas=local_maximas, mask_shape=mask_shape
        )
        if ind == 0:
            all_jet_core_mask = current_jet_core_mask
            continue
        all_jet_core_mask = np.dstack([all_jet_core_mask, current_jet_core_mask])
    return all_jet_core_mask


def create_mask_using_local_maximas(local_maximas, mask_shape):
    """
    Will create a mask with the same dimensions as the inputted mask_shape
    and only the local maxima as value 1. All other values will be 0.

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    local_maximas : np.array
        An array containing values relating to the index of the local maximas in a mask with shape: mask_shape

    mask_shape : tuple
        Values to create a mask of given shape from.

    Returns
    ----------
    mask : np.array
        An array with shape: mask_shape with only the indexes of the local maximas as 1 All other values will 0.

    """
    empty_mask = np.zeros(mask_shape)
    for i, j in local_maximas:
        empty_mask[i, j] = 1
    return empty_mask


def get_all_jet_region_contour_mask(row, local_maximas_dict):
    """
    Runs 'get_jet_region_contour_mask' on each longitude in data with one time step

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    row : xr.DataArray
        Data of single time unit containing the variables: 'potential_jet_cores', and the coordinates: 'lon', 'lat', 'plev'

    local_maximas : np.array
        An array containing values relating to the index of the local maximas in a mask with shape: mask_shape

    Returns
    ----------
    all_jet_region_contour_mask : np.array
        A 3-D array the indexes (lat, plev) of the local maximas for each longitude as 1 All other values will 0.

    """
    all_jet_region_contour_mask = np.array([])
    for ind, lon in enumerate(row.lon):
        current_jet_region_contour_mask = get_jet_region_contour_mask(
            row["potential_jet_regions"].sel(lon=lon), local_maximas_dict[float(lon)]
        )
        if ind == 0:
            all_jet_region_contour_mask = current_jet_region_contour_mask
            continue
        all_jet_region_contour_mask = np.dstack(
            [all_jet_region_contour_mask, current_jet_region_contour_mask]
        )
    return all_jet_region_contour_mask


def get_jet_region_contour_mask(potential_jet_regions, local_maximas):
    """
    Will create a mask based on regions (above ws_threshold) around the jet cores that are
    contiguous.
    The mask will contain categorical values (e.g. 1, 2) for each cluster of jet regions
    contain a jet core. All other values will be 0 (i.e not a jet core or jet region)

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    potential_jet_regions : xr.DataArray
        Windspeeds of regions above the jet region threshold

    local_maximas : np.array
        An array containing values relating to the index of the local maximas (should be the same shape as input: potential_jet_regions)

    Returns
    ----------
    jet_region_contour_mask : np.array
        A 2-D array with shape: mask_shape with only the indexes of the local maximas as 1 All other values will 0.

    """
    potential_jet_regions_mask, _ = scipy.ndimage.label(potential_jet_regions)
    actual_jet_region_nums = get_jet_region_numbers(
        local_maximas, potential_jet_regions_mask
    )
    jet_region_contour_mask = subset_jet_region_mask_to_regions_with_cores(
        potential_jet_regions_mask, actual_jet_region_nums
    )

    # Refine jet regions to new values (i.e. if labels have been removed this resets them)
    jet_region_contour_mask, _ = scipy.ndimage.label(jet_region_contour_mask)
    return jet_region_contour_mask


def get_jet_region_numbers(local_maximas, potential_jet_regions_mask):
    """
    Will only return clusters ID numbers of jet regions with an actual jet core in them.

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    local_maximas : np.array
        An array containing values relating to the index of the local maximas (should be the same shape as input: potential_jet_regions)

    potential_jet_regions_mask : numpy.array
        Array of jet regions clusters as returned by scipy.ndimage.label

    Returns
    ----------
    actual_jet_region_nums : list
        A list of cluster ID numbers with local maximas in them

    """
    actual_jet_region_nums = []
    for i, j in local_maximas:
        valid_region = potential_jet_regions_mask[i, j]
        if valid_region not in actual_jet_region_nums:
            actual_jet_region_nums.append(valid_region)
    return actual_jet_region_nums


def subset_jet_region_mask_to_regions_with_cores(
    potential_jet_regions_mask, actual_jet_region_nums
):
    """
    Will subset jet region mask to only those with an actual jet core in them.

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    potential_jet_regions_mask : numpy.array
        Array of jet regions clusters as returned by scipy.ndimage.label

    actual_jet_region_nums : list
        A list of cluster ID numbers with local maximas in them

    Returns
    ----------
    potential_jet_regions_mask : numpy.array
        Subset version of jet regions mask containing only jet region clusters with local maxima in them

    """
    for reg_num in np.unique(potential_jet_regions_mask):
        if reg_num not in actual_jet_region_nums:
            potential_jet_regions_mask[potential_jet_regions_mask == reg_num] = 0
    return potential_jet_regions_mask


def get_all_jet_regions_mask(row, all_jet_region_contour_mask, local_maximas_dict):
    """
    Get all jet regions mask by looping over each longitude in row.
    See also docs of 'refine_jet_region_to_leftright_and_abovebelow'

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    row : xr.DataArray
        Data of single time unit containing the coordinates: 'lon'

    all_jet_region_contour_mask : np.array
        A 3-D array the indexes (lat, plev) of the local maximas for each longitude as 1 All other values will 0.

    local_maximas_dict : dict
        Keys are longitude coordinates, values are local maxima locations

    Returns
    ----------
    all_jet_regions_mask : np.array
        A refined version of the inputted 'all_jet_region_contour_mask' of the local maximas for each longitude as 1 All other values will 0.

    """
    all_jet_regions_mask = np.array([])  # initialise jet region array
    for ind, lon in enumerate(row.lon):
        current_jet_region_contour_mask_one_lon = all_jet_region_contour_mask[
            ::, ::, ind
        ]
        current_jet_region_mask = np.zeros_like(
            current_jet_region_contour_mask_one_lon
        )  # create empty mask

        current_local_maximas = local_maximas_dict[float(lon)]
        for local_maxima in current_local_maximas:
            refined_result = refine_jet_region_to_leftright_and_abovebelow(
                current_jet_region_contour_mask_one_lon,
                local_maxima[0],
                local_maxima[1],
            )
            current_jet_region_mask += refined_result
        if ind == 0:
            all_jet_regions_mask = current_jet_region_mask
            continue
        all_jet_regions_mask = np.dstack(
            [all_jet_regions_mask, current_jet_region_mask]
        )
    return all_jet_regions_mask


def refine_jet_region_to_leftright_and_abovebelow(array, x, y):
    """
    This method will remove all values not left, right, above or below the input
    x, y coordinates in a 2-D array.

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    array : np.array
        A 2-D array with numeric values relating to the jet region (i.e. dtypes float or int)
    x : int
        x-coordinate of a local maxima to subset jet region by
    y: int
        y-coordinate of a local maxima to subset input array by

    Returns
    ----------
    refined_array : np.array
        New 2-D array with values only around x and y coordinates

    """
    rows, cols = array.shape
    refined_array = np.zeros_like(array)

    if 0 <= x < rows and 0 <= y < cols and array[x, y] != 0:
        # Refine the current cluster
        refined_array[x, y] = array[x, y]
        left_col = y - 1
        while left_col >= 0 and array[x, left_col] != 0:
            refined_array[x, left_col] = array[x, left_col]
            left_col -= 1

        right_col = y + 1
        while right_col < cols and array[x, right_col] != 0:
            refined_array[x, right_col] = array[x, right_col]
            right_col += 1

        up_row = x - 1
        while up_row >= 0 and array[up_row, y] != 0:
            refined_array[up_row, y] = array[up_row, y]
            up_row -= 1

        down_row = x + 1
        while down_row < rows and array[down_row, y] != 0:
            refined_array[down_row, y] = array[down_row, y]
            down_row += 1

    return refined_array


def get_values_along_a_line_between_two_coordinates(data, start_point, end_point):
    """
    Get all values along a shortest path between two coordinates in a 2-D numpy array.

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    arr : np.array
        A 2-D array with numeric values (i.e. dtypes float or int)

    Returns
    ----------
    local_maxima : np.array
        A 2-D array with values relating to the index of the local maximas

    Examples
    ----------
    .. code-block:: python

        import numpy as np

        # Example array
        data = np.array([
            [1, 2, 3, 4, 5],
            [2, 3, 4, 5, 4],
            [3, 4, 5, 6, 7],
            [4, 5, 6, 7, 6],
            [5, 6, 7, 8, 9]
        ])

        local_maxima_indices = find_local_maxima(data)

        for i, j in local_maxima_indices:
            print(f"Local maximum at ({i}, {j}): {data[i, j]}")

    """
    # Calculate the differences in coordinates
    dx = end_point[1] - start_point[1]
    dy = end_point[0] - start_point[0]

    # Calculate the number of steps required for the line
    num_steps = max(abs(dx), abs(dy))

    # Calculate the step sizes for each coordinate
    x_step = dx / num_steps
    y_step = dy / num_steps

    # Initialize lists to store the coordinates along the path
    path_coordinates = []
    for step in range(num_steps + 1):
        row = int(round(start_point[0] + step * y_step))
        col = int(round(start_point[1] + step * x_step))
        path_coordinates.append((row, col))

    # Extract values along the path
    values_along_path = [float(data[row, col]) for row, col in path_coordinates]

    return values_along_path


def has_ws_drop_between_cores(ws_between_cores, ws_drop_threshold):
    """
    Will check for a windspeed drop of a given threshold between cores

    Parameters
    ----------
    ws_between_cores : np.array
        A 1-D array with windspeed between two cores as determined by 'get_values_along_a_line_between_two_coordinates'
    ws_drop_threshold : int or float
        Windspeed threshold to check

    Returns
    ----------
    out : boolean
        True if windspeeds input have a drop greater than 'ws_drop_threshold'

    """
    if any(ws_between_cores[0] - np.array(ws_between_cores) > ws_drop_threshold):
        return True
    elif any(ws_between_cores[-1] - np.array(ws_between_cores) > ws_drop_threshold):
        return True
    else:
        return False


def get_index_of_cores_to_drop_based_on_multicore_regions(
    local_maxima_ind_slices, current_local_maximas, multi_core_region_ws
):
    """
    Get index of cores to remove in same region based on whether they contain the first occurence of
    the maximum windspeed or not.

    Parameters
    ----------
    local_maxima_ind_slices : np.array of np.arrays
        Collection of index slices for current local maximas
    current_local_maximas : np.array
        Collection of indexes of current local maxima
    multi_core_region_ws : xr.DataArray
        Data containing region windspeeds

    Returns
    ----------
    index_of_cores_to_drop : np.array
        Mask of jet cores in current longitude plev/lat slice with multi core regions formatted
    """
    index_of_cores_to_drop = []
    # Loop through each local_maxima_slice
    for local_maxima_ind_slice in local_maxima_ind_slices:
        if len(local_maxima_ind_slice) <= 1:
            # This would indicate a single core, so will not be removed
            continue
        # Get current max windspeed for region (as fallback). Will select first occurence of max.
        index_of_ws_maxima = get_current_region_ws_maxima_lat_and_plev_ind(
            current_local_maximas, local_maxima_ind_slice, multi_core_region_ws
        )
        previous_local_maxima = [0, 0]  # temporary
        for ind, local_maxima_ind in enumerate(local_maxima_ind_slice):
            current_local_maxima = current_local_maximas[int(local_maxima_ind)]
            if ind == 0:
                previous_local_maxima = current_local_maxima
                continue
            # Check if current or previous is the maximum (so will need to be retained)
            current_is_index_of_max = np.array_equal(
                np.array([current_local_maxima[0], current_local_maxima[1]]),
                index_of_ws_maxima,
            )
            previous_is_index_of_max = np.array_equal(
                np.array([previous_local_maxima[0], previous_local_maxima[1]]),
                index_of_ws_maxima,
            )
            # Store index of cores to drop
            if current_is_index_of_max:
                index_of_cores_to_drop.append(
                    [previous_local_maxima[0], previous_local_maxima[1]]
                )
            else:
                index_of_cores_to_drop.append(
                    [current_local_maxima[0], current_local_maxima[1]]
                )
                if not previous_is_index_of_max:
                    # occurs only if there are 3+ cores in region and maximum is not between current and previous
                    index_of_cores_to_drop.append(
                        [previous_local_maxima[0], previous_local_maxima[1]]
                    )
            previous_local_maxima = current_local_maxima
    return index_of_cores_to_drop


def get_current_region_ws_maxima_lat_and_plev_ind(
    current_local_maximas, local_maxima_ind_slice, multi_core_region_ws
):
    """
    Get the plev and lat index of current region windspeed maxima.
    Will select the first occurence of maxima, if there is more than one (rare).

    Parameters
    ----------
    current_local_maximas : np.array
        Collection of indexes of current local maxima
    local_maxima_ind_slice : np.array
        Index slice for current local maxima
    multi_core_region_ws : xr.DataArray
        Data containing region windspeeds

    Returns
    ----------
    output : np.array
        Array of the plev and latitude
    """
    all_current_local_maximas_ind = current_local_maximas[
        list(local_maxima_ind_slice.astype(int))
    ]
    all_current_local_maximas_data = multi_core_region_ws.isel(
        plev=all_current_local_maximas_ind[:, 0],
        lat=all_current_local_maximas_ind[:, 1],
    )
    if all_current_local_maximas_data.isnull().all():
        return []
    max_ws_indices = np.unravel_index(
        np.nanargmax(all_current_local_maximas_data),
        all_current_local_maximas_data.shape,
    )
    current_max_data = all_current_local_maximas_data[
        max_ws_indices[0], max_ws_indices[1]
    ]
    # select first occurence of maximum if multiple
    if current_max_data.size > 1:
        current_max_data = current_max_data[0]
    iplev = list(multi_core_region_ws.plev.values).index(current_max_data.plev)
    ilat = list(multi_core_region_ws.lat.values).index(current_max_data.lat)
    return np.array([iplev, ilat])


def run_checks_on_jet_cores_and_return_jet_cores(
    row,
    initial_jet_core_masks,
    local_maximas_dict,
    jet_core_lat_distance,
    ws_drop_threshold,
):
    """
    This method runs two checks on the jet cores to check whether there are regions with multiple jet cores.
    Firstly, it checks whether regions with multiple jet cores are more than a certain distance apart
    (default is 15 degrees, see 'jet_core_lat_distance'), and hence seperate cores.
    Secondly, it will check whether the windspeed between two cores drops below a threshold
    (default is 25 m/s, see 'ws_drop_threshold'), if so it will remove the latter core.

    Component of method of algorithm originally introduced in Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    row : xarray.Dataset
        Data of single time unit containing the variables: 'jet_region_contour_mask', and the coordinates: 'lon', 'lat', 'plev'
    initial_jet_core_masks :
        Initial mask of jet cores to check.
    ws_drop_threshold : int or float
        Threshold for drop in windspeed along the line between cores (default: 25 m/s)
    jet_core_lat_distance : int or float
        Threshold for maximum distance between cores to be counted the same (default: 15 degrees)

    Returns
    ----------
    jet_core_masks : numpy.array
        Final jet cores mask of shape of 'initial_jet_core_masks'.

    """
    jet_core_masks = np.copy(initial_jet_core_masks)
    for lon_ind, lon in enumerate(row.lon):
        jet_region_contour_one_lon = row["jet_region_contour_mask"].sel(lon=lon)
        ws_one_lon = row["ws"].sel(lon=lon)
        current_local_maximas = local_maximas_dict[float(lon)]
        # Get first occurence of the maximum windspeed
        index_of_first_maximum = np.dstack(np.where(ws_one_lon == ws_one_lon.max()))[0]

        core_and_location = []
        for core_ind, local_maxima in enumerate(current_local_maximas):
            region_within = jet_region_contour_one_lon[
                local_maxima[0], local_maxima[1]
            ]  # this is the region contour that the jet core is found within
            core_and_location.append([core_ind, float(region_within)])
        core_and_location = np.array(core_and_location)

        if len(core_and_location) == 0:
            # No cores found, so selecting the cell with maximum wind speed for this lon slice
            core_and_location = index_of_first_maximum

        region_ind, num_cores_in_region = np.unique(
            core_and_location[::, 1], return_counts=True
        )
        multi_core_regions = region_ind[num_cores_in_region > 1]

        if len(multi_core_regions) == 0:
            # no multi-core regions found
            continue

        for multi_core_region in multi_core_regions:
            # get all the local maxima within regions of multi cores
            local_maxima_inds = core_and_location[::, 0][
                core_and_location[::, 1] == multi_core_region
            ]

            # Check 1. Test if cores are more than 15 degrees away
            previous_lat = None
            for local_maxima_ind in local_maxima_inds:
                local_maxima = current_local_maximas[int(local_maxima_ind)]
                if not previous_lat:
                    # set previous latitude to check for latitude distance in next
                    previous_lat = local_maxima[1]
                else:
                    current_lat = local_maxima[1]
                    if (
                        abs(row["lat"][previous_lat] - row["lat"][current_lat])
                        > jet_core_lat_distance
                    ):
                        previous_lat = current_lat
                        continue

            # Check 2. Test if cores have ws drop between them
            multi_core_region_ws = ws_one_lon.where(
                jet_region_contour_one_lon == multi_core_region
            )
            for multi_core_region in multi_core_regions:
                # get all the local maxima within regions of multi cores
                local_maxima_inds = core_and_location[::, 0][
                    core_and_location[::, 1] == multi_core_region
                ]
                # Step 1 of 2: loop over local maxima in multicore region and find if they need to split further
                slicepoints_in_region = []
                previous_local_maxima = [0, 0]  # temporary
                for ind, local_maxima_ind in enumerate(local_maxima_inds):
                    current_local_maxima = current_local_maximas[int(local_maxima_ind)]
                    if ind == 0:
                        previous_local_maxima = current_local_maxima
                        continue
                    windspeeds_between_cores = (
                        get_values_along_a_line_between_two_coordinates(
                            multi_core_region_ws,
                            start_point=previous_local_maxima,
                            end_point=current_local_maxima,
                        )
                    )
                    # Get conditions for dropping cores
                    ws_drops_below_threshold = has_ws_drop_between_cores(
                        windspeeds_between_cores, ws_drop_threshold=ws_drop_threshold
                    )
                    if ws_drops_below_threshold:
                        slicepoints_in_region.append(ind)
                    previous_local_maxima = current_local_maxima
                # Step 2 of 2: Remove cores in same region and retain core with maximum windspeed
                local_maxima_ind_slices = data_utils.slice_array_by_index_breaks(
                    local_maxima_inds, slicepoints_in_region
                )
                index_of_cores_to_drop = (
                    get_index_of_cores_to_drop_based_on_multicore_regions(
                        local_maxima_ind_slices,
                        current_local_maximas,
                        multi_core_region_ws,
                    )
                )
                for index_of_core_to_drop in index_of_cores_to_drop:
                    jet_core_masks[
                        index_of_core_to_drop[0], index_of_core_to_drop[1], lon_ind
                    ] = 0
    return jet_core_masks


def get_empty_local_wind_maxima_data(
    data, expected_dims=("time", "plev", "lat", "lon")
):
    """
    Will add a new data var of zeros for local wind maxima

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    data : xarray.Dataset
        Data which should containing the variables: 'ua' and 'va', and the expected dims e.g.: 'lon', 'lat', 'plev' and 'time'.

    expected_dims : tuple
        Expected dimensions of input data (default: "time", "plev", "lat", "lon")

    Returns
    ----------
    data : xarray.Dataset
        Data containing zeros array of (time, plev, lat, lon) dimensions
    """
    expected_dims_set = set(expected_dims)
    data_dim_set = set(data.dims)
    actual_dims = tuple(expected_dims_set.intersection(data_dim_set))
    dim_sizes = tuple(len(data[coord]) for coord in actual_dims)
    data["local_wind_maxima"] = (
        actual_dims,
        np.zeros(dim_sizes),
    )
    return data


def get_potential_local_wind_maximas_by_ws_threshold(ws_slice, ws_threshold=30):
    """
    Will return a 2-d array of potential local windspeed maximas

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    ws_slice : xarray.Dataset
        Data slice of windspeed that has only 'lat' and 'lon' dims

    ws_threshold : int or float
        windspeed threshold to apply (default=30 ms-1)
    Returns
    ----------
    ws_slice : xarray.Dataset
        Data slice of windspeed (lat, lon only) with ws_threshold applied
    """
    return ws_slice.where(lambda x: x > ws_threshold).fillna(0.0)


def get_local_wind_maxima_by_timeunit(row, ws_threshold):
    """
    Get local wind maxima by timeunit (i.e. day)

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305 who originally use 30 m/s as their ws threshold

    Parameters
    ----------
    row : xarray.Dataset
        Data of single time unit containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev'

    ws_threshold : int or float
        windspeed threshold to apply

    Returns
    ----------
    row : xarray.Dataset
        Data of a single time unit containing 0 or 1 value (local wind maxima) for that time unit
    """
    if "local_wind_maxima" not in row.data_vars:
        raise ValueError("local_wind_maxima needs to be defined.")

    row = row.transpose("plev", "lat", ...)
    for lon in row["lon"]:
        current = row.sel(lon=lon, method="nearest")
        pot_local_maximas = get_potential_local_wind_maximas_by_ws_threshold(
            current["ws"], ws_threshold
        ).data
        ind_local_wind_maximas = data_utils.get_local_maxima(pot_local_maximas, axis=1)
        # Turn into 2-d numpy array
        ind_local_wind_maximas = np.array(
            [
                [arr1, arr2]
                for arr1, arr2 in zip(
                    ind_local_wind_maximas[0], ind_local_wind_maximas[1]
                )
            ]
        )
        for plev_ind, lat_ind in ind_local_wind_maximas:
            row["local_wind_maxima"].loc[
                dict(
                    plev=current["plev"].data[plev_ind],
                    lon=lon,
                    lat=current["lat"].data[lat_ind],
                )
            ] = 1.0
    return row


def get_number_of_timeunits_per_monthyear_with_local_wind_maxima(data):
    """
    Will resample by each month and return number of timeunits (i.e. day) with local wind maxima.

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    data : xarray.Dataset
        Input data with local wind maxima values by time unit (i.e. day)

    Returns
    ----------
    data : xarray.Dataset
        Data containing zeros array of (time, plev, lat, lon) dimensions

    """
    data["local_wind_maxima_by_monthyear"] = (
        data["local_wind_maxima"]
        .resample(time="MS")
        .sum()
        .rename({"time": "monthyear"})
    )
    return data


def subdivide_local_wind_maxima_into_stj_pfj(
    data, local_wind_column_name="local_wind_maxima_by_monthyear"
):
    """
    Subdivide the local_wind_maxima values into the Subtropical Jet (STJ) and Polar Front Jet (PFJ) based on Table 1 pg. 2709
    from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305.
    After the method in that paper, categorisation for the Northern Hemisphere STJ is only possible in DJF,
    and the latitude bands used are not based on December, March, June or September in Table 1 of that study.

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    data : xarray.Dataset
        Input data with local wind maxima values

    Returns
    ----------
    data : xarray.Dataset
        data with polar_front_jet and subtropical_jet subdivisions
    """
    DJF_STJ = data.sel(
        monthyear=data.monthyear.dt.month.isin([12, 1, 2]), lat=slice(15, 40)
    )[local_wind_column_name]
    MAM_SON_PFJ = data.sel(
        monthyear=data.monthyear.dt.month.isin([3, 4, 5, 9, 10, 11]),
        lat=slice(10, 70),
    )[local_wind_column_name]
    JJA_PFJ = data.sel(
        monthyear=data.monthyear.dt.month.isin([6, 7, 8]), lat=slice(30, 60)
    )[local_wind_column_name]
    SH_STJ = data.sel(lat=slice(-40, -15))[local_wind_column_name]
    SH_PFJ = data.sel(lat=slice(-70, -41))[local_wind_column_name]
    data["polar_front_jet"] = xr.merge([MAM_SON_PFJ, JJA_PFJ, SH_PFJ])[
        local_wind_column_name
    ]
    data["subtropical_jet"] = xr.merge([DJF_STJ, SH_STJ])[local_wind_column_name]
    return data


def run_jet_occurence_and_centre_alg_on_one_day(row, occurence_ws_threshold):
    """
    Runs JetStreamCoreIdentificationAlgorithm method on a single day.

    Component of method from Kuang et al (2014) https://doi.org/10.1007/s00704-013-0994-x

    Parameters
    ----------
    row : xarray.Dataset
        Data of single time unit containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev'
    occurence_ws_threshold : int or float
        Threshold used to identify a jet-stream occurence point

    Returns
    ----------
    occ_alg.output_data : xarray.Dataset
        Data with jet occurence and centre points (1 for occurence, 2 for centre)
    """
    # Step 0. Squeeze time for method
    if "time" in row.dims:
        row = row.squeeze("time")

    # Step 1. Get jet occurences
    row["jet_occurence"] = row["ws"].where(
        lambda val: (val >= occurence_ws_threshold), 0
    )

    # Step 2. Get jet centers (local maximas) at each longitude
    local_maximas_dict = get_local_maximas_at_each_longitude(
        row, check_diagonals=True, var_name="jet_occurence"
    )

    # Step 3. Loop through the local maximas and make an mask of jet centers
    mask_shape = row["ws"].isel(lon=0).shape
    all_jet_centers_mask = get_all_jet_core_mask(
        local_maximas_dict=local_maximas_dict, mask_shape=mask_shape
    )

    # Step 4. Clip to between 0-1 and return row
    row["jet_occurence"] = np.clip(row["jet_occurence"], 0, 1)

    row["jet_centers"] = (
        ("plev", "lat", "lon"),
        np.clip(all_jet_centers_mask, 0, 1),
    )
    return row


class JetStreamOccurenceAndCentreAlgorithm:
    """
    Jet-stream occurence and centre algorithm.

    Component of method from Kuang et al (2014) https://doi.org/10.1007/s00704-013-0994-x
    """

    def __init__(self, data, occurence_ws_threshold=30):
        """
        Parameters
        ----------
        data : xarray.Dataset
           Data of single time unit containing the variables: 'ua' and 'va', and the only the coordinates: 'lon', 'lat'
        occurence_ws_threshold : int or float
            Threshold used to identify a jet-stream occurence point
        """
        try:
            assert occurence_ws_threshold > 0
        except Exception as e:
            raise ValueError(
                "Occurence wind-speed threshold needs to be more than 0"
            ) from e

        # Load in data as a pressure level 2d wind-speed slice
        self.plev_ws_slice = windspeed_utils.PressureLevelWindSpeedSlice(data).values
        self.plev_ws_slice["jet_ocurrence1_jet_centre2"] = self.plev_ws_slice[
            "ws"
        ].copy()
        self.plev_ws_slice["jet_ocurrence1_jet_centre2"] = self.plev_ws_slice[
            "jet_ocurrence1_jet_centre2"
        ].where(lambda x: x >= occurence_ws_threshold)
        self._jet_occurence = self.plev_ws_slice
        self._lat_resolution = abs(
            float(self.plev_ws_slice["lat"][1] - self.plev_ws_slice["lat"][0])
        )
        self._lon_resolution = abs(
            float(self.plev_ws_slice["lon"][1] - self.plev_ws_slice["lon"][0])
        )

        # make_duplicate data for output
        self.output_data = self._jet_occurence.copy(deep=True)

        # Initialise lists needed for search algorithm
        self._all_coords = []
        self._lats_with_3 = []
        self._lats_for_search = []
        self._jet_centres = []

        # Needed to keep track of algorithm
        self.algorithm_has_run = False

    @classmethod
    def run_algorithm(cls, data):
        """
        Class method for running the algorithm.
        """
        alg_output = cls(data)
        alg_output.run()
        return alg_output

    def run(self):
        """
        Run the algorithm.
        """
        self._get_all_coords_of_jet_occurence()
        self._all_coords_arr = np.array(self._all_coords)
        self._calc_jet_centre_points()
        self._get_jet_centre_data()
        self._label_jet_occurence()

        self.algorithm_has_run = True

    def _get_jet_centre_data(self):
        """
        Calculates jet-stream centres based on if one jet-stream occurence grid
        is surrounded by 8 cells of jet-stream occurence (default is 30 m/s)
        """
        # to fix: there's got to be a quicker way
        for centre in self._jet_centres:
            self.output_data["jet_ocurrence1_jet_centre2"].loc[
                dict(lat=centre[0], lon=centre[1])
            ] = 2

    def _get_all_coords_of_jet_occurence(self):
        """
        Get all coords with the jet occurence value.
        """
        for val in self._jet_occurence["jet_ocurrence1_jet_centre2"].notnull():
            if val.any():
                for sub_val in val:
                    if sub_val:
                        self._all_coords.append(
                            [float(sub_val["lat"]), float(sub_val["lon"])]
                        )

    def _get_all_latitudes_available_in_3by3_grid(self):
        """
        Get all latitudes available in 3 by 3 grid.
        """
        for lat in self._lats_with_3:
            if (
                lat - self._lat_resolution in self._lats_with_3
                and lat + self._lat_resolution in self._lats_with_3
            ):
                self._lats_for_search.append(lat)

    def _calc_jet_centre_points(self):
        """
        Will return a list of the coordinates for all jet-centre points
        """
        for coord in self._all_coords_arr:
            coord_ws = float(
                self.output_data.sel(lat=coord[0], lon=coord[1], method="nearest")["ws"]
            )
            lat_grid_vals = np.arange(
                coord[0] - self._lat_resolution,
                coord[0] + self._lat_resolution + 0.01,
                self._lat_resolution,
            )
            lat_grid_vals = lat_grid_vals[
                (lat_grid_vals >= float(self.output_data["lat"].min()))
                & (lat_grid_vals <= float(self.output_data["lat"].max()))
            ]
            lon_grid_vals = np.arange(
                coord[1] - self._lon_resolution,
                coord[1] + self._lon_resolution + 0.01,
                self._lon_resolution,
            )
            lon_grid_vals = lon_grid_vals % 360  # loop lon around
            matrix_vals_to_check = np.array(
                np.meshgrid(lat_grid_vals, lon_grid_vals)
            ).T.reshape(-1, 2)
            add_coord = True
            for val_to_check in matrix_vals_to_check:
                if val_to_check[-1] not in self.output_data["lon"]:
                    continue
                if (
                    float(
                        self.output_data.sel(
                            lat=val_to_check[0], lon=val_to_check[1], method="nearest"
                        )["ws"]
                    )
                    > coord_ws
                ):
                    add_coord = False
                    break
            if add_coord:
                self._jet_centres.append(coord)

    def _label_jet_occurence(self):
        """
        Will label all non 2 values of windspeed for the occurence
        Used in Kuang et al. 2014
        """
        self.output_data["jet_ocurrence1_jet_centre2"] = self.output_data[
            "jet_ocurrence1_jet_centre2"
        ].where(lambda x: np.isfinite(x), 0)
        self.output_data["jet_ocurrence1_jet_centre2"] = self.output_data[
            "jet_ocurrence1_jet_centre2"
        ].where(lambda x: ((x == 0) | (x == 2)), 1)


def run_jet_core_algorithm_on_one_day(row, ws_core_threshold, ws_boundary_threshold):
    """
    Runs JetStreamCoreIdentificationAlgorithm method on a single time unit.

    Component of method of jet_core_identification_algorithm in jsmetrics and is inspired by
    the method from Manney et al. (2011) (https://doi.org/10.5194/acp-11-6115-2011),
    which is described in Section 3.1 of that study.

    Parameters
    ----------
    row : xarray.Dataset
        Data of single time unit containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev'

    ws_core_threshold : int or float
        Threshold used for jet-stream core point
    ws_boundary_threshold : int or float
        Threshold for jet-stream boundary point

    Returns
    ----------
    row : xarray.Dataset
        Data for one time unit containing jet-cores (ID number relates to each unique core)
    """
    row["jet_core_id"] = (
        ("plev", "lat", "lon"),
        np.zeros((row["plev"].size, row["lat"].size, row["lon"].size)),
    )
    for lon in row["lon"]:
        current = row.sel(lon=lon, method="nearest")
        core_alg = JetStreamCoreIdentificationAlgorithm(
            current,
            ws_core_threshold=ws_core_threshold,
            ws_boundary_threshold=ws_boundary_threshold,
        )
        core_alg.run()
        row["jet_core_id"].loc[dict(lon=lon)] = core_alg.output_data["core_id"]
    return row


class JetStreamCoreIdentificationAlgorithm:
    """
    Jet-stream core identification algorithm.

    Component of method of jet_core_identification_algorithm in jsmetrics and is inspired by
    the method from Manney et al. (2011) (https://doi.org/10.5194/acp-11-6115-2011),
    which is described in Section 3.1 of that study.

    Methods
    -------
    run:
        run algorithm
    run_algorithm:
        class method for running algorithm
    """

    def __init__(self, data, ws_core_threshold=40, ws_boundary_threshold=30):
        """
        Input will need to be longitudinal slice of windspeed values.

        Component of method of jet_core_identification_algorithm in jsmetrics and is inspired by
        the method from Manney et al. (2011) (https://doi.org/10.5194/acp-11-6115-2011),
        which is described in Section 3.1 of that study.


        Parameters
        ----------
        data : xarray.Dataset
            Data of single time unit containing the variables: 'ua' and 'va', and the coordinates: 'lon', 'lat', 'plev'
        ws_core_threshold : int or float
            Threshold used for jet-stream core point (default=40)
        ws_boundary_threshold : int or float
            Threshold for jet-stream boundary point (default=30)
        """
        try:
            assert (
                ws_core_threshold > ws_boundary_threshold
                and ws_core_threshold > 0
                and ws_boundary_threshold > 0
            )
        except Exception as e:
            raise ValueError(
                "Windspeed core threshold needs to be more than boundary\
                    threshold and both need to be more than 0"
            ) from e
        # Transpose data
        data = data.transpose(*(..., "lat", "plev"))

        # Step 1. make windspeed slice
        self._lat_ws_slice = windspeed_utils.LatitudeWindSpeedSlice(data)

        # Step 2. Get core and potential boundary points
        self._labelled_data = self._lat_ws_slice.label_slice(
            self._lat_ws_slice["ws"] < ws_core_threshold, "Core"
        )
        self._labelled_data = self._labelled_data.where(
            (self._lat_ws_slice["ws"] < ws_boundary_threshold)
            | (self._lat_ws_slice["ws"] > ws_core_threshold),
            other="Potential Boundary",
        )

        # Step 3. Get indexes of jet-stream cores and potential boundaries
        (
            self._initial_core_ids,
            self._pot_boundary_ids,
        ) = self._get_indexes_of_core_and_boundaries()

        # Step 4. Set variables needed to keep track of
        self._current_core_lat = -1
        self._currently_a_core = None
        self.algorithm_has_run = False

    def __repr__(self):
        """
        Representation of the class. Have it return the labelled data
        """
        if not self.algorithm_has_run:
            print(
                "A total of %d initial Jet-stream cores have been found\
                 in the wind-speed slice"
                % (self._labelled_data["ws"].where(lambda x: x == "Core").count())
            )
            print(
                "A total of %d potential Jet-stream boundaries have been found\
                 in the wind-speed slice"
                % (
                    self._labelled_data["ws"]
                    .where(lambda x: x == "Potential Boundary")
                    .count()
                )
            )
            return repr(self._labelled_data)
        else:
            print("A total of %d cores have been discovered" % (self.num_of_cores))
            return repr(self.output_data)

    @classmethod
    def run_algorithm(cls, data, ws_core_threshold=40, ws_boundary_threshold=30):
        """
        Class method for running algorithm
        """
        js_algorithm = cls(
            data,
            ws_core_threshold=ws_core_threshold,
            ws_boundary_threshold=ws_boundary_threshold,
        )

        js_algorithm.run()
        return js_algorithm

    def run(self):
        """
        Runs algorithm.
        """
        self.final_jet_cores = self._get_jet_core_boundary()
        self.output_data = self._add_cores_to_data()
        self.algorithm_has_run = True

    def _get_indexes_of_core_and_boundaries(self):
        """
        Will return the indexes in the ws data that ARE jet-stream cores
        and COULD BE jet-stream core boundaries
        """
        pot_boundary_ids = np.where(self._labelled_data["ws"] == "Potential Boundary")
        initial_core_ids = np.where(self._labelled_data["ws"] == "Core")
        pot_boundary_ids = np.stack(pot_boundary_ids, axis=-1)
        initial_core_ids = np.stack(initial_core_ids, axis=-1)
        return initial_core_ids, pot_boundary_ids

    @staticmethod
    def _get_indexes_to_check(pot_boundary):
        """
        Will return an array of indexes to check for potential boundaries
        or jetstream cores.
        """
        vals_to_check = []
        if pot_boundary[0] != 0:
            vals_to_check.append([pot_boundary[0] - 1, pot_boundary[1]])
        vals_to_check.append([pot_boundary[0] + 1, pot_boundary[1]])
        if pot_boundary[1] != 0:
            vals_to_check.append([pot_boundary[0], pot_boundary[1] - 1])
        vals_to_check.append([pot_boundary[0], pot_boundary[1] + 1])
        return vals_to_check

    def _get_pot_jetcore_area(self, vals, area, core_found=False):
        """
        Recursive function that will return the IDs of a jet core boundary
        i.e. above 30 m/s surrounding a core of 40 m/s
        Will check one an area of potential boundaries contains a core
        and thus can be called boundaries.
        """
        vals_copy = vals.copy()
        for val in vals:
            if val in area:
                continue
            if val in self._initial_core_ids.tolist():
                core_found = True
                # look for a new core if it is more than 15 degrees away
                if (
                    val[1] - self._current_core_lat > 15
                    and val[1] > self._current_core_lat
                ):
                    #  print('THIS IS A NEW CORE')
                    self._current_core_lat = val[1]
                area.append(val)
                new_vals = self._get_indexes_to_check(val)
                vals_copy.extend(new_vals)
                vals_copy = data_utils.remove_duplicates(vals_copy)
                return self._get_pot_jetcore_area(
                    vals_copy, area=area, core_found=core_found
                )

            elif val in self._pot_boundary_ids.tolist():
                area.append(val)
                new_vals = self._get_indexes_to_check(val)
                vals_copy.extend(new_vals)
                vals_copy = data_utils.remove_duplicates(vals_copy)
                return self._get_pot_jetcore_area(
                    vals_copy, area=area, core_found=core_found
                )
            else:
                vals_copy.remove(val)
                continue

        # reset current core variables
        self._current_core_lat = -1
        self._currently_a_core = None
        return area, core_found

    def _get_jet_core_boundary(self):
        """
        Recursive function that will return the IDs of all jet core boundaries
        i.e. above 30 m/s surrounding a core of 40 m/s
        Will check if an area of potential boundaries contains a core
        and thus can be called boundaries.
        """
        already_covered = []
        js_core_indexes = []
        id_number = 0
        for pot_boundary in self._pot_boundary_ids:
            if pot_boundary.tolist() in already_covered:
                continue
            vals_to_check = self._get_indexes_to_check(pot_boundary)
            area, core_found = self._get_pot_jetcore_area(vals_to_check, area=[])
            already_covered.extend(area)
            already_covered = data_utils.remove_duplicates(already_covered)
            # add area to js_core_indexes if part of core
            if core_found:
                id_number += 1
                js_core_indexes.extend([{"id": id_number, "index_of_area": area}])
        self.num_of_cores = id_number
        return js_core_indexes

    def _add_cores_to_data(self):
        """
        Add cores to data.
        """
        self._lat_ws_slice.values["core_id"] = (
            ("plev", "lat"),
            np.zeros(
                (
                    self._lat_ws_slice.values["plev"].size,
                    self._lat_ws_slice.values["lat"].size,
                )
            ),
        )
        for jet_core in self.final_jet_cores:
            for lat, plev in jet_core["index_of_area"]:
                self._lat_ws_slice.values["core_id"].loc[
                    dict(
                        lat=self._lat_ws_slice.values["lat"].data[lat],
                        plev=self._lat_ws_slice.values["plev"].data[plev],
                    )
                ] = jet_core["id"]
        return self._lat_ws_slice.values
