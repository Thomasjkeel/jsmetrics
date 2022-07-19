# -*- coding: utf-8 -*-

"""
    Algorithms and calculations for the metrics used to
    identify or classify jet-stream in the literature.

    This file is in order of publish year of the metrics
    (see details_for_all_metrics.py)
"""

# imports
# import collections
import numpy as np

# import matplotlib.pyplot
# import xarray as xr
# import scipy.fftpack
# import scipy.interpolate
# # import shapely.geometry
from . import data_utils, windspeed_utils

# docs
__author__ = "Thomas Keel"
__email__ = "thomas.keel.18@ucl.ac.uk"
__status__ = "Development"


def get_sum_weighted_ws(data, all_plevs_hPa):
    """
    Get sum of weighted windspeed.
    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    sum weighted windspeed = integral(p2, p1)(u^2+v^2)^(1/2)dp
    where p1, p2 is min, max pressure level

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
        raise TypeError(
            "array of pressure level needs to be list or numpy.array"
        )

    sum_weighted_ws = 0
    for plev, (i, plev_hPa) in zip(data["plev"], enumerate(all_plevs_hPa)):
        if i != 0:
            plev_hPa = plev_hPa - all_plevs_hPa[i - 1]
        sum_weighted_ws += (
            (data.sel(plev=plev)["ua"] ** 2 + data.sel(plev=plev)["va"] ** 2)
            ** (1 / 2)
        ) * plev_hPa
    return sum_weighted_ws


def get_weighted_average_ws(sum_weighted_ws, all_plevs_hPa):
    """
    Component of method from Koch et al (2006) https://doi.org/10.1002/joc.1255

    weighted average windspeed = 1/(p2-p1) * sum average windspeed
    where p1, p2 is min, max pressure level

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
    if not isinstance(all_plevs_hPa, (list, np.ndarray)):
        raise TypeError(
            "array of pressure level needs to be a list or numpy.array"
        )
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
    ):
        raise ValueError("Plev units need to be mbar, Pa or hPa")

    plevs = np.array([plev for plev in data["plev"]])
    if data["plev"].units == "Pa":
        plevs = plevs / 100
    return plevs


def get_local_jet_maximas_by_timeunit_by_plev(row):
    """
    Component of method from Schiemann et al 2009 https://doi.org/10.1175/2008JCLI2625.1
    TODO: add checks
    TODO: will only work is 1 day is the resolution
    TODO: maybe combine with pena-ortiz method

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
    )  # TODO
    for lon in row["lon"]:
        for plev in row["plev"]:
            current = row.sel(lon=lon, plev=plev)
            current = current.where(
                (abs(current["ws"]) >= 30) & (current["ua"] > 0)
            )
            local_maxima_lat_inds = data_utils.get_local_maxima(
                current["ws"].data
            )[0]
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


def run_jet_core_algorithm_on_one_day(
    row, ws_core_threshold, ws_boundary_threshold
):
    """
    Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
    Runs JetStreamCoreIdentificationAlgorithm method on a single time unit

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
    )  # TODO
    for lon in row["lon"]:
        current = row.sel(lon=lon)
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
    Jet-stream core identification algorithm
    Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
    """

    def __init__(self, data, ws_core_threshold=40, ws_boundary_threshold=30):
        """
        Component of method  from Manney et al. (2011) https://doi.org/10.5194/acp-11-6115-2011
        input will need to be longitudinal slice of windspeed values


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

    def __add__(self, other):
        print("TODO: need to implement behaviour for adding slices together")
        pass

    def __repr__(self):
        """
        Representation of the class. Have it return the labelled data
        """
        if not self.algorithm_has_run:
            print(
                "A total of %d initial Jet-stream cores have been found\
                 in the wind-speed slice"
                % (
                    self._labelled_data["ws"]
                    .where(lambda x: x == "Core")
                    .count()
                )
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
            print(
                "A total of %d cores have been discovered"
                % (self.num_of_cores)
            )
            return repr(self.output_data)

    @classmethod
    def run_algorithm(
        cls, data, ws_core_threshold=40, ws_boundary_threshold=30
    ):
        """
        Class method for running algorithm
        """
        js_algorithm = cls(
            data,
            ws_core_threshold=ws_core_threshold,
            ws_boundary_threshold=ws_boundary_threshold,
        )

        js_core_indexes = js_algorithm.run()
        return js_core_indexes

    def run(self):
        self.final_jet_cores = self._get_jet_core_boundary()
        self.output_data = self._add_cores_to_data()
        self.algorithm_has_run = True

    def _get_indexes_of_core_and_boundaries(self):
        """
        Will return the indexes in the ws data that ARE jet-stream cores
        and COULD BE jet-stream core boundaries
        """
        pot_boundary_ids = np.where(
            self._labelled_data["ws"] == "Potential Boundary"
        )
        initial_core_ids = np.where(self._labelled_data["ws"] == "Core")
        pot_boundary_ids = np.stack(pot_boundary_ids, axis=-1)
        initial_core_ids = np.stack(initial_core_ids, axis=-1)
        return initial_core_ids, pot_boundary_ids

    @staticmethod
    def _get_indexes_to_check(pot_boundary):
        """
        Will return an array of indexes to check for potential boundaries
        or jetstream cores.
        Used in Manney et al. 2011
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
            area, core_found = self._get_pot_jetcore_area(
                vals_to_check, area=[]
            )
            already_covered.extend(area)
            already_covered = data_utils.remove_duplicates(already_covered)
            # add area to js_core_indexes if part of core
            if core_found:
                id_number += 1
                js_core_indexes.extend(
                    [{"id": id_number, "index_of_area": area}]
                )
        self.num_of_cores = id_number
        return js_core_indexes

    def _add_cores_to_data(self):
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
