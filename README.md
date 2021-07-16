# jet-stream-metrics

To-Do:
- way of dealing with data from different sources (some sort of data translator module or maybe included in tests)
    - for example what if 'v' or 'v-wind' is passed to func instead of 'va'
    - for example what if 'mbar' or 'model levels' instead of 'plev'
- add raise exception statements to required 
- Add way of subsetting longitude if it wraps around 0-360
- Improve readme(s) and add links (DOI?)m (https://github.com/othneildrew/Best-README-Template/blob/master/README.md)
- add lat and lon subsetting to jetstream_metrics_dict
- look into bash arrays and running multiple metrics
- add benchmark speed tests for each metric based on standard sizes of the dataset
- argparse to parse a list for additional arguments (do more user testing)
- add doctests to classes?


#### Why class based algorithms:
- see https://softwareengineering.stackexchange.com/questions/294076/how-to-represent-an-algorithm-as-a-class
