from . import jetstream_metrics

# RULES for JETSTREAM_METRICS dictionary:
    # 1. must have the keys: 'variables', 'coords' and 'metric'
    # 2. 'variables' will contain the required CMIP6? model output variable names
    # 3. 'coords' will contain the required CMIP6? standard coords and each coord will provide a list of 2 values: mininum value for coord and maximum value for coord and must be a number
    #    3.1 for 'plev' coord it is in millibars and higher pressure is considered the minimum value (e.g. 85000 mbar - 50000 mbar) 
    # 4. 'metric' is the name of a function in jetstream_metrics.py

JETSTREAM_METRIC_DICT = {
    "Koch2006":
        {"variables": ["ua","va"], "coords": {"plev": [10000, 40000]}, "metric": jetstream_metrics.koch_et_al_2006, "description":"Koch et al. 2006"},
    "Woolings2010":
         {"variables": ["ua"], "coords": {"plev": [70000, 92500]}, "metric": jetstream_metrics.woolings_et_al_2010, "description":"Woolings et al. 2010 TODO"}, # , "exact_coords": {"plev": [92500, 85000, 77500, 70000]}
    "ScreenSimmonds2013":
        {"variables":["zg"], "coords":{"plev": [50000, 50000]}, "metric": jetstream_metrics.screen_and_simmonds_2013, "description":"Screen & Simmonds 2013 TODO"},
    "Kuang2014":
        {"variables":["ua", "va"], "coords":{"plev": [20000,25000]}, "metric": jetstream_metrics.kuang_et_al_2014, "description":"Kuang et al. 2014. NOTE: Should be 200hPa only"},
    "FrancisVavrus2015":
        {"variables":["ua","va"], "coords":{"plev": [50000,50000]}, "metric": jetstream_metrics.francis_vavrus_2015, "description":"Francis & Vavrus 2015"},
    "LocalWaveActivty":
        {"variables":["zg"], "coords":{"plev": [50000,50000]}, "metric": jetstream_metrics.local_wave_activity, "description":"First introduced by Huang & Nakamura 2015, then modified by Chen et al. 2015 and Martineau et al. 2017"},
    "Cattiaux2016":
        {"variables":["zg"], "coords":{"plev": [50000, 50000]}, "metric": jetstream_metrics.cattiaux_et_al_2016, "description":"Cattiaux et al. 2016"},
    "Ceppi2018":
        {"variables":["ua"], "coords":{"plev": [85000, 85000], "lat": [30,60]}, "metric": jetstream_metrics.ceppi_et_al_2018, "description":"Ceppi et al. 2018"},
    "Kern2018":
        {"variables":["ua", "va"], "coords":{"plev": [0, 100000]}, "metric": jetstream_metrics.kern_et_al_2018, "description":"Kern 2018"},
    "Simpson2018":
        {"variables":["zg"], "coords":{"plev": [70000, 70000]}, "metric": jetstream_metrics.simpson_et_al_2018, "description":"Simpson et al. 2018"},
    "ChemkeMing2020":
        {"variables":["ua","va"], "coords":{"plev": [0, 100000]}, "metric": jetstream_metrics.chemke_and_ming_2020, "description":"Chemke & Ming 2020"}
    } 

