import argparse
from experiments import test_uvwind

ALL_EXPERIMENTS = {'experiment1':test_uvwind.main}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run an experiment.')
    parser.add_argument('--experiment', type=str, help='path to a python script')
    args = parser.parse_args()
    print("Running script: ", args.experiment)

    ALL_EXPERIMENTS[args.experiment]()
