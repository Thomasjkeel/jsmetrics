from jsmetrics.metrics import jet_core_algorithms, jet_statistics, waviness_metrics


# RULES for METRIC_DETAILS dictionary:
# 1. must have the keys: 'variables', 'coords' 'metric', 'plev_units', 'metric', 'name', 'description, and 'doi'
# 2. 'variables' will contain the required variable names in the data conforming to the CMIP Controlled Vocabulary (Taylor et al. 2011)
# 3. 'coords' will contain the required conforming to the CMIP Controlled Vocabulary
#       and each coord will provide a list of 2 integer values: a minimum and maximum value
#    3.1 for 'plev' coord it is in Pa/mbar and higher pressure is considered
#        the minimum value (e.g. 85000 Pmbar - 50000 mbar)
# 4. 'plev_units' is units for the pressure level (i.e. Pa, mbar, hPa)
# 5. 'metric' is the relative path to the function
# 6. 'name' is the name of the paper
# 7. 'doi' is the DOI link


METRIC_DETAILS = {
    "Koch2006": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 40000]},
        "plev_units": "Pa",
        "metric": jet_core_algorithms.koch_et_al_2006,
        "name": "Koch et al. 2006",
        "description": """This a two-step method for detecting \'jet-event occurences\'.
                         based on a calculating a weighted average windspeed and
                        then applying a windspeed threshold to isolate jet events from
                        the weighted average. Note that the original method also includes
                        a third step (production of a climatology), but we instead provide
                        an example of how to calculate this third step in the 'Examples'
                        section of the methods docstring. For more information see
                        'Notes' for this method available on the ReadTheDocs for jsmetrics""",
        "doi": "https://doi.org/10.1002/joc.1255",
    },
    "ArcherCaldeira2008": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 40000]},
        "plev_units": "Pa",
        "metric": jet_statistics.archer_caldeira_2008,
        "name": "Archer & Caldeira 2008",
        "description": "This method extracts three monthly-averaged jet stream properties\
                        via integrated quantities (windspeed, pressure and latitude) from u-component wind.\
                        For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics",
        "doi": "https://doi.org/10.1029/2008GL033614",
    },
    "Schiemann2009": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 50000]},
        "plev_units": "Pa",
        "metric": jet_core_algorithms.schiemann_et_al_2009,
        "name": "Schiemann et al 2009",
        "description": "This method detects jet occurrences, whereby each jet occurence is\
                        detected based on three rules applied to inputted wind speed. These\
                        rules are based on whether input windspeed is a local maxima\
                        (in y and z dimension), and applying various thresholds.\
                        Originally the threshold is >= 30 ms for wind vector and \
                        >= 0 ms for input u-wind. This method is described in section 2\
                        of that study. For more information see 'Notes' for this\
                        method available on the ReadTheDocs for jsmetrics",
        "doi": "https://doi.org/10.1175/2008JCLI2625.1",
    },
    "Woollings2010_NorthAtlantic": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 92500], "lat": [15, 75], "lon": [120, 180]},
        "plev_units": "Pa",
        "metric": jet_statistics.woollings_et_al_2010,
        "name": "Woollings et al. 2010 North Atlantic",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics",
        "doi": "https://onlinelibrary.wiley.com/doi/10.1002/qj.625",
    },
    "Manney2011": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 100000]},
        "plev_units": "Pa",
        "metric": jet_core_algorithms.manney_et_al_2011,
        "name": "Manney et al. 2011",
        "description": "This method detects jet cores and defines a boundary region\
                        beside those cores based on two windspeed thresholds.\
                        Two additional checks are applied after initial detection\
                        of cores to check whether cores within the same windspeed region\
                        are part of the same feature (see docstring for this method for\
                        more information).",
        "doi": "https://doi.org/10.5194/acp-11-6115-2011",
    },
    "BarnesPolvani2013_NorthAtlantic": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 85000], "lat": [0, 90], "lon": [300, 360]},
        "plev_units": "Pa",
        "metric": jet_statistics.barnes_polvani_2013,
        "name": "Barnes & Polvani 2013 North Atlantic",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1175/JCLI-D-12-00536.1",
    },
    "BarnesPolvani2013_NorthPacific": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 85000], "lat": [0, 90], "lon": [135, 235]},
        "plev_units": "Pa",
        "metric": jet_statistics.barnes_polvani_2013,
        "name": "Barnes & Polvani 2013 North Pacific",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1175/JCLI-D-12-00536.1",
    },
    "BarnesPolvani2013_SouthernHemisphere": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 85000], "lat": [-90, 0]},
        "plev_units": "Pa",
        "metric": jet_statistics.barnes_polvani_2013,
        "name": "Barnes & Polvani 2013 Southern Hemisphere",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1175/JCLI-D-12-00536.1",
    },
    "PenaOrtiz2013": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 40000]},
        "plev_units": "Pa",
        "metric": jet_core_algorithms.penaortiz_et_al_2013,
        "name": "Pena-Ortiz et al. 2013",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1002/jgrd.50305",
    },
    # "ScreenSimmonds2013_NorthAmerica_NorthAtlantic": {
    #     "variables": ["zg"],
    #     "coords": {"plev": [50000, 50000], "lon": [220, 360]},
    #     "plev_units": "Pa",
    #     "metric": jetstream_metrics.screen_and_simmonds_2013,
    #     "name": "Screen & Simmonds 2013 North America and North Atlantic",
    #     "description": "",
    #     "doi": "https://doi.org/10.1002/grl.50174",
    # },
    # "ScreenSimmonds2013_Europe_Atlantic": {
    #     "variables": ["zg"],
    #     "coords": {"plev": [50000, 50000], "lon": [300, 420/60]},
    #     "plev_units": "Pa",
    #     "metric": jetstream_metrics.screen_and_simmonds_2013,
    #     "name": "Screen & Simmonds 2013 Europe and Atlantic",
    #     "description": "",
    #     "doi": "https://doi.org/10.1002/grl.50174",
    # },
    # "ScreenSimmonds2013_NorthAmerica_WestPacific": {
    #     "variables": ["zg"],
    #     "coords": {"plev": [50000, 50000], "lon": [180, 300]},
    #     "plev_units": "Pa",
    #     "metric": jetstream_metrics.screen_and_simmonds_2013,
    #     "name": "Screen & Simmonds 2013 North America and West Pacific",
    #     "description": "",
    #     "doi": "https://doi.org/10.1002/grl.50174",
    # },
    # "ScreenSimmonds2013_Asia_EastPacific": {
    #     "variables": ["zg"],
    #     "coords": {"plev": [50000, 50000], "lon": [60, 180]},
    #     "plev_units": "Pa",
    #     "metric": jetstream_metrics.screen_and_simmonds_2013,
    #     "name": "Screen & Simmonds 2013 Asia and East Pacific",
    #     "description": "",
    #     "doi": "https://doi.org/10.1002/grl.50174",
    # },
    "Kuang2014": {
        "variables": ["ua", "va"],
        "coords": {"plev": [20000, 20000]},
        "plev_units": "Pa",
        "metric": jet_core_algorithms.kuang_et_al_2014,
        "name": "Kuang et al. 2014.",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1007/s00704-013-0994-x",
    },
    "BarnesPolvani2015_NorthAmerica_NorthAtlantic": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 92500], "lat": [30, 70], "lon": [230, 350]},
        "plev_units": "Pa",
        "metric": jet_statistics.barnes_polvani_2015,
        "name": "Barnes & Polvani 2015",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1175/JCLI-D-14-00589.1",
    },
    "FrancisVavrus2015": {
        "variables": ["ua", "va"],
        "coords": {"plev": [50000, 50000]},
        "plev_units": "Pa",
        "metric": waviness_metrics.francis_vavrus_2015,
        "name": "Francis & Vavrus 2015",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1088/1748-9326/10/1/014005",
    },
    "Cattiaux2016": {
        "variables": ["zg"],
        "coords": {"plev": [50000, 50000]},
        "plev_units": "Pa",
        "metric": waviness_metrics.cattiaux_et_al_2016,
        "name": "Cattiaux et al. 2016",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1002/2016GL070309",
    },
    "GrisePolvani2016": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [-65, -30]},
        "plev_units": "Pa",
        "metric": jet_statistics.grise_polvani_2016,
        "name": "Grise & Polvani 2016",
        "description": "This method calculates the latitude of the midlatitude eddy-driven jet\
                        by finding the peak value of the input u-wind field.\
                        A polynomial fit is then applied to get an appropriate value of 'jet_lat'\
                        at a resolution 0.01 degrees. For more information see 'Notes' for this\
                        method available on the ReadTheDocs for jsmetrics.",
        "doi": "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1002/2015JD024687",
    },
    "BarnesSimpson2017_NorthAtlantic": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 70000], "lat": [1, 90], "lon": [280, 350]},
        "plev_units": "Pa",
        "metric": jet_statistics.barnes_simpson_2017,
        "name": "Barnes & Simpson 2017 North Atlantic",
        "description": "This method defines two outputs: 'jet_lat' and 'jet_speed' \
                        which are defined as the latitude and speed of the 10-day-averaged\
                        maximum zonally-averaged wind speed. For more information see 'Notes'\
                        for this method available on the ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1175/JCLI-D-17-0299.1",
    },
    "BarnesSimpson2017_NorthPacific": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 70000], "lat": [1, 90], "lon": [120, 230]},
        "plev_units": "Pa",
        "metric": jet_statistics.barnes_simpson_2017,
        "name": "Barnes & Simpson 2017 North Pacific",
        "description": "This method defines two outputs: 'jet_lat' and 'jet_speed' \
                        which are defined as the latitude and speed of the 10-day-averaged\
                        maximum zonally-averaged wind speed. For more information see 'Notes'\
                        for this method available on the ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1175/JCLI-D-17-0299.1",
    },
    "Bracegirdle2018_SouthernHemisphere": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [-75, -10]},
        "plev_units": "Pa",
        "metric": jet_statistics.bracegirdle_et_al_2018,
        "name": " Bracegirdle et al. 2018",
        "description": "This method calculates the seasonal and annual jet-stream position ('JPOS')\
                        and strength ('JSTR') by applying a 0.075 degree cubic spline interpolation to zonally-averaged\
                        wind climatology and selecting the maximum. For more information see 'Notes'\
                        for this method available on the ReadTheDocs for jsmetrics.",
        "doi": "https://doi.org/10.1175/JCLI-D-17-0320.1",
    },
    "Ceppi2018_NorthAtlantic_Europe": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [30, 60], "lon": [300, 60]},
        "plev_units": "Pa",
        "metric": jet_statistics.ceppi_et_al_2018,
        "name": "Ceppi et al. 2018 North Atlantic",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "10.1175/JCLI-D-17-0323.1",
    },
    "Ceppi2018_NorthPacific": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [30, 60], "lon": [140, 240]},
        "plev_units": "Pa",
        "metric": jet_statistics.ceppi_et_al_2018,
        "name": "Ceppi et al. 2018 North Pacific",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "10.1175/JCLI-D-17-0323.1",
    },
    "Ceppi2018_SouthernHemisphere": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [-60, -30]},
        "plev_units": "Pa",
        "metric": jet_statistics.ceppi_et_al_2018,
        "name": "Ceppi et al. 2018 Southern Hemisphere",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "10.1175/JCLI-D-17-0323.1",
    },
    "Zappa2018_NorthAtlantic": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [20, 70], "lon": [300, 360]},
        "plev_units": "Pa",
        "metric": jet_statistics.zappa_et_al_2018,
        "name": "Zappa et al. 2018 North Atlantic",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "10.1002/2017GL076096",
    },
    "Zappa2018_NorthPacific": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [20, 70], "lon": [140, 240]},
        "plev_units": "Pa",
        "metric": jet_statistics.zappa_et_al_2018,
        "name": "Zappa et al. 2018 North Pacific",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "10.1002/2017GL076096",
    },
    "Kerr2020_NorthernHemisphere": {
        "variables": ["ua"],
        "coords": {"plev": [50000, 50000], "lat": [20, 70]},
        "plev_units": "Pa",
        "metric": jet_statistics.kerr_et_al_2020,
        "name": "Kerr et al. 2020",
        "description": "For more information see 'Notes' for this method available on the\
                        ReadTheDocs for jsmetrics.",
        "doi": "10.1029/2020JD032735",
    },
    # "BlackportFyfe2022_NorthAtlantic": {
    #     "variables": ["ua"],
    #     "coords": {"plev": [70000, 70000], "lat": [15, 75], "lon": [300, 360]},
    #     "plev_units": "Pa",
    #     "metric": jet_statistics.blackport_fyfe_2022,
    #     "name": "Blackport & Fyfe 2022",
    #     "description": "",
    #     "doi": "10.1126/sciadv.abn3112",
    # },
}
