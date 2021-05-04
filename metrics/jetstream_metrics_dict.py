from . import jetstream_metrics

# RULES for JETSTREAM_METRICS dictionary:
    # 1. must have the keys: 'variables', 'coords' and 'metric'
    # 2. 'variables' will contain the required CMIP6? model output variable names
    # 3. 'coords' will contain the required CMIP6? standard coords and each coord will provide a list of 2 values: mininum value for coord and maximum value for coord and must be a number
    #    3.1 for 'plev' coord it is in millibars and higher pressure is considered the minimum value (e.g. 85000 mbar - 50000 mbar) 
    # 4. 'metric' is the name of a function in jetstream_metrics.py

JETSTREAM_METRIC_DICT = {
    "Koch2006":
        {"variables": ["ua","va"], "coords": {"plev": [40000, 10000]}, "metric": jetstream_metrics.koch_et_al_2006, "description":"Koch et al. 2006"},
    "Woolings2010":
         {"variables": ["ua"], "coords": {"plev": [92500,  70000]}, "metric": jetstream_metrics.woolings_et_al_2010, "description":"Woolings et al. 2010 TODO"}, # , "exact_coords": {"plev": [92500, 85000, 77500, 70000]}
    "Hudson2012":
        {"variables":["ozone"], "coords":{"DU": [450, 250]}, "metric": jetstream_metrics.hudson_2012, "description":"Hudson et al. 2012"},
    "ScreenSimmonds2014":
        {"variables":["zg"], "coords":{"plev": [50000, 50000]}, "metric": jetstream_metrics.screen_and_simmonds_2014, "description":"Screen & Simmonds 2014"},
    "Kuang2014":
        {"variables":["ua"], "coords":{"plev": [20000,20000]}, "metric": jetstream_metrics.kuang_et_al_2014, "description":"Kuang et al. 2014"},
    "HuangNakamura2015":
        {"variables":["rv850?"], "coords":{"plev": [85000,85000]}, "metric": jetstream_metrics.huang_and_nakamura_2015, "description":"Huang & Nakamura 2015"},
    "Cattiaux2016":
        {"variables":["zg"], "coords":{"plev": [50000, 50000]}, "metric": jetstream_metrics.cattiaux_et_al_2016, "description":"Cattiaux et al. 2016"},
    "Ceppi2018":
        {"variables":["ua"], "coords":{"plev": [85000, 85000]}, "metric": jetstream_metrics.ceppi_et_al_2018, "description":"Ceppi et al. 2018"},
    "Kern2018":
        {"variables":["ua", "va"], "coords":{"plev": [100000,0]}, "metric": jetstream_metrics.kern_et_al_2018, "description":"Kern 2018"},
    "Simpson2018":
        {"variables":["zg"], "coords":{"plev": [70000, 70000]}, "metric": jetstream_metrics.simpson_et_al_2018, "description":"Simpson et al. 2018"},
    "ChemkeMing2020":
        {"variables":["ua","va"], "coords":{"plev": [100000,0]}, "metric": jetstream_metrics.chemke_and_ming_2020, "description":"Chemke & Ming 2020"}
    } 
