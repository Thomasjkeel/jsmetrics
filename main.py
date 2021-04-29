import argparse
from experiments import test_uvwind

ALL_EXPERIMENTS = {'experiment1':test_uvwind.main}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run an experiment.')
    parser.add_argument('-e','--experiment', type=str, help='path to a python script')
    parser.add_argument('-d', '--data', type=str, help='path to a netcdf4 data')
    args = parser.parse_args()
    print("Running script: ", args.experiment)

    ## TODO: make error handling smarter
    try:
        ALL_EXPERIMENTS[args.experiment]()
    except Exception as e:
        print(e, 'is not a known experiment')
