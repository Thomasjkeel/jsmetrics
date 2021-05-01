import argparse
from experiments import test_available_metrics

ALL_EXPERIMENTS = {'test_metrics':test_available_metrics.main}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run an experiment.')
    parser.add_argument('-e','--experiment', type=str, help='path to a python script')
    parser.add_argument('-d', '--data', type=str, help='path to a netcdf4 data')
    args = parser.parse_args()
    print("Running script: ", args.experiment)

    ## TODO: make error handling smarter
    if args.experiment not in ALL_EXPERIMENTS.keys():
        print(args.experiment, 'is not a known experiment')
    else:
        ALL_EXPERIMENTS[args.experiment]()
        
