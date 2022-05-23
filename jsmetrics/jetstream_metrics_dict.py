from . import jetstream_metrics


# RULES for JETSTREAM_METRICS dictionary:
# 1. must have the keys: 'variables', 'coords' 'metric' and 'name'
# 2. 'variables' will contain the required CMIP6? model output variable names
# 3. 'coords' will contain the required CMIP6? standard coords
#       and each coord will provide a list of 2 values: mininum value
#       for coord and maximum value for coord and must be a number
#    3.1 for 'plev' coord it is in Pa/mbar and higher pressure is considered
#        the minimum value (e.g. 85000 Pmbar - 50000 mbar)
# 4. 'plev_units' is units for the pressure level
# 5. 'metric' is the name of a function in jetstream_metrics.py
# 6. 'name' is the paper name
# 7. 'doi' is the DOI link


JETSTREAM_METRIC_DICT = {
    "Koch2006": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 40000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.koch_et_al_2006,
        "name": "Koch et al. 2006",
        "description": "",
        "doi": "https://doi.org/10.1002/joc.1255",
    },
    "ArcherCaldeira2008": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 40000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.archer_caldeira_2008,
        "name": "Archer & Caldeira 2008",
        "description": "",
        "doi": "https://doi.org/10.1029/2008GL033614",
    },
    # "BartonEllis2009": {
    #     "variables": ["ua", "va"],
    #     "coords": {"plev": [30000, 30000]},
    #     "plev_units": "Pa",
    #     "metric": jetstream_metrics.barton_ellis_2009,
    #     "name": "Barton & Ellis 2009",
    #     "description": "",
    #     "doi": "https://doi.org/10.1002/joc.1750",
    # },
    "Schiemann2009": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 50000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.schiemann_et_al_2009,
        "name": "Schiemann et al 2009",
        "description": "",
        "doi": "https://doi.org/10.1175/2008JCLI2625.1",
    },
    "Woolings2010": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 92500], "lat": [15, 75], "lon": [120, 180]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.woolings_et_al_2010,
        "name": "Woolings et al. 2010",
        "description": "",
        "doi": "https://onlinelibrary.wiley.com/doi/10.1002/qj.625",
    },  # , "exact_coords": {"plev": [92500, 85000, 77500, 70000]}
    "Manney2011": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 40000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.manney_et_al_2011,
        "name": "Manney et al. 2011",
        "description": "surface to 0.1 hPa",
        "doi": "https://doi.org/10.5194/acp-11-6115-2011",
    },
    "BarnesPolvani2013": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 85000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.barnes_polvani_2013,
        "name": "Barnes & Polvani 2013",
        "description": "",
        "doi": "https://doi.org/10.1175/JCLI-D-12-00536.1",
    },
    "PenaOrtiz2013": {
        "variables": ["ua", "va"],
        "coords": {"plev": [10000, 40000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.penaortiz_et_al_2013,
        "name": "Pena-Ortiz et al. 2013",
        "description": "",
        "doi": "https://doi.org/10.1002/jgrd.50305",
    },
    "ScreenSimmonds2013": {
        "variables": ["zg"],
        "coords": {"plev": [50000, 50000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.screen_and_simmonds_2013,
        "name": "Screen & Simmonds 2013",
        "description": "",
        "doi": "https://doi.org/10.1002/grl.50174",
    },
    "Kuang2014": {
        "variables": ["ua", "va"],
        "coords": {"plev": [20000, 25000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.kuang_et_al_2014,
        "name": "Kuang et al. 2014.",
        "description": "NOTE: Should be 200hPa only",
        "doi": "https://doi.org/10.1007/s00704-013-0994-x",
    },
    "BarnesPolvani2015": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 92500], "lat": [30, 70], "lon": [230, 350]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.barnes_polvani_2015,
        "name": "Barnes & Polvani 2015",
        "description": "",
        "doi": "https://doi.org/10.1175/JCLI-D-14-00589.1",
    },
    "FrancisVavrus2015": {
        "variables": ["ua", "va"],
        "coords": {"plev": [50000, 50000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.francis_vavrus_2015,
        "name": "Francis & Vavrus 2015",
        "description": "",
        "doi": "https://doi.org/10.1088/1748-9326/10/1/014005",
    },
    # "LocalWaveActivity":
    #     {"variables": ["zg"], "coords":{"plev": [50000, 50000]},
    #       "metric": jetstream_metrics.local_wave_activity,
    #       "description": "First introduced by Huang & Nakamura 2015, then
    #                 modified by Chen et al. 2015 and Martineau et al. 2017"},
    "Cattiaux2016": {
        "variables": ["zg"],
        "coords": {"plev": [50000, 50000]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.cattiaux_et_al_2016,
        "name": "Cattiaux et al. 2016",
        "description": "",
        "doi": "https://doi.org/10.1002/2016GL070309",
    },
    "BarnesSimpson2017": {
        "variables": ["ua"],
        "coords": {"plev": [70000, 70000], "lat": [0, 90], "lon": [110, 350]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.barnes_simpson_2017,
        "name": "Barnes & Simpson 2017",
        "description": "",
        "doi": "https://doi.org/10.1175/JCLI-D-17-0299.1",
    },
    "GrisePolvani2017": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [-65, -30]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.grise_polvani_2017,
        "name": "Grise & Polvani 2017",
        "description": "",
        "doi": "https://journals.ametsoc.org/doi/10.1175/JCLI-D-16-0849.1",
    },
    # "Molnos2017":
    #     {"variables": ["ua", "va"], "coords":{"plev":[50000, 15000]},
    #      "metric": jetstream_metrics.molnos_et_al_2017,
    #      "description": "Molnos et al 2017"},
    "Bracegirdle2018": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [-75, -10]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.bracegirdle_et_al_2018,
        "name": " Bracegirdle et al. 2018",
        "description": "",
        "doi": "https://doi.org/10.1175/JCLI-D-17-0320.1",
    },
    "Ceppi2018": {
        "variables": ["ua"],
        "coords": {"plev": [85000, 85000], "lat": [30, 60]},
        "plev_units": "Pa",
        "metric": jetstream_metrics.ceppi_et_al_2018,
        "name": "Ceppi et al. 2018",
        "description": "",
        "doi": "10.1175/JCLI-D-17-0323.1",
    },
    # "Kern2018":
    #     {"variables": ["ua", "va"], "coords":{"plev": [0, 100000]}, "plev_units": "Pa",
    #      "metric": jetstream_metrics.kern_et_al_2018,
    #      "description": "Kern 2018"},
    # "Rikus2018":
    #     {"variables": ["ua"], "coords":{"plev": [0, 100000]}, "plev_units": "Pa",
    #      "metric": jetstream_metrics.rikus_2018,
    #      "description": "Rikus 2018"},
    # "ChemkeMing2020":
    #     {"variables": ["ua","va"], "coords":{"plev": [0, 100000]}, "plev_units": "Pa",
    #      "metric": jetstream_metrics.chemke_and_ming_2020,
    #      "description": "Chemke & Ming 2020"}
}
