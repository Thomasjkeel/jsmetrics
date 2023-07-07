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
import xarray as xr

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_sum_weighted_ws(data, all_plevs_hPa):
    """
    Get sum of weighted windspeed.
    sum weighted windspeed = integral(p2, p1)(u^2+v^2)^(1/2)dp
    where p1, p2 is min, max pressure level

    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    Parameters
    ----------
    data : xarray.Dataset
        Data containing u- and v-component wind
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
    """
    Calculates weighted average wind-speed:
    weighted average windspeed = 1/(p2-p1) * sum average windspeed
    where p1, p2 is min, max pressure level

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

    plevs = np.array([plev for plev in data["plev"]])
    if data["plev"].units == "Pa":
        plevs = plevs / 100
    return plevs


def get_local_jet_maximas_by_oneday_by_plev(row):
    """
    Get local jet maxima for one day.

    Component of method from Schiemann et al 2009 https://doi.org/10.1175/2008JCLI2625.1
    NOTE: will only work if 1 day is the resolution

    Parameters
    ----------
    row : xarray.Dataset
        Data of a sinlge time unit containing windspeed (ws), plev, lat, lon

    Returns
    ----------
    row : xarray.Dataset
        Data of a sinlge time unit with value for jet-maxima (1 == maxima, 0 == none)

    """
    row["jet_maxima"] = (
        ("plev", "lat", "lon"),
        np.zeros((row["plev"].size, row["lat"].size, row["lon"].size)),
    )
    for lon in row["lon"]:
        for plev in row["plev"]:
            current = row.sel(lon=lon, plev=plev, method="nearest")
            current = current.where((abs(current["ws"]) >= 30) & (current["ua"] > 0))
            local_maxima_lat_inds = data_utils.get_local_maxima(current["ws"].data)[0]
            if len(local_maxima_lat_inds) > 0:
                for lat_ind in local_maxima_lat_inds:
                    row["jet_maxima"].loc[
                        dict(
                            lat=current["lat"].data[lat_ind],
                            lon=lon,
                            plev=plev,
                        )
                    ] = 1.0
    return row


def run_jet_core_algorithm_on_one_day(row, ws_core_threshold, ws_boundary_threshold):
    """
    Runs JetStreamCoreIdentificationAlgorithm method on a single time unit.

    Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

    Parameters
    ----------
    row : xarray.Dataset
        Data of single time unit containing u- and v-component wind
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

    Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

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

        Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011

        Parameters
        ----------
        data : xarray.Dataset
            Data of single time unit containing u- and v-component wind
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


def get_empty_local_wind_maxima_data(
    data, expected_dims=("time", "plev", "lat", "lon")
):
    """
    Will add a new data var of zeros for local wind maxima

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    data : xarray.Dataset
        Input data with (time, plev, lat, lon) dimensions

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
        Data slice of windspeed that has only lat and lon dims

    ws_threshold : int or float
        windspeed threshold to apply (default=30 ms-1)
    Returns
    ----------
    ws_slice : xarray.Dataset
        Data slice of windspeed (lat, lon only) with ws_threshold applied
    """
    return ws_slice.where(lambda x: x > ws_threshold).fillna(0.0)


def get_local_wind_maxima_by_timeunit(row):
    """
    Get local wind maxima by timeunit (i.e. day)

    Component of method from Pena-Ortiz (2013) https://doi.org/10.1002/jgrd.50305

    Parameters
    ----------
    row : xarray.Dataset
        Data of a single time unit containing windspeed (ws)

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
            current["ws"], 30
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
    Subdivide the local_wind_maxima values into the Subtropical Jet (STJ) and Polar Front Jet (PFJ) based on pg. 2709.

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
        monthyear=data.monthyear.dt.month.isin([1, 2, 12]), lat=slice(15, 40)
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
        Single time-unit of dataset containing u- and v-component wind
    occurence_ws_threshold : int or float
        Threshold used to identify a jet-stream occurence point

    Returns
    ----------
    occ_alg.output_data : xarray.Dataset
        Data with jet occurence and centre points (1 for occurence, 2 for centre)
    """
    occ_alg = JetStreamOccurenceAndCentreAlgorithm(row, occurence_ws_threshold)
    occ_alg.run()
    return occ_alg.output_data


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
            Single time-unit of dataset containing u- and v-component wind. Needs to be 2-dimensions (lat, lon)
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
        # TODO: there's got to be a quicker way
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
