import argparse
from experiments import get_available_metrics, calculate_metrics_from_data
from tests import test_plot

ALL_EXPERIMENTS = {'get_available_metrics':{"description": "Will print out all metrics available for this dataset","script":get_available_metrics.main},
                   'calc_metrics':{"description":"will calcualte any number of metrics. Use -m tag to state which one", "script":calculate_metrics_from_data.main},
                   'test_plot': {"description": "Will test the subset and calculate a metric from the data and produce a plot (under experiments/figures)", "script": test_plot.main}}



def list_all_experiments():
    """
        Prints out all the experiments and their shorthand names
        TODO: maybe right allign this print statement
    """
    print("\nAll experiments:\n")
    print("-"*20)
    for exp in ALL_EXPERIMENTS.keys():
        print(exp, "—", ALL_EXPERIMENTS[exp]["description"])
    print("-"*20)
    print("\n")
    return


def run_experiments(args):
    """
        Will run experiment off of users -e (experiment tag) input
        TODO: make error handling smarter
    """
    if not args['experiment']:
        print('Error: No experiment declared. Use: \'-e experiment_name\' or \'-l to list all experiments\'')
    elif args['experiment'] not in ALL_EXPERIMENTS.keys():
        print('Error: \'%s\' is not a known experiment. Run \'-l\' for list of experiments' % (args['experiment']))
    else:
        print("Running script: ", args['experiment'])
        if args['metrics']:
            print("Metrics to use: ", args['metrics'])
            # TODO: add allowance for zero to multiple metrics
            ALL_EXPERIMENTS[args['experiment']]["script"](**args)
        else:
            ALL_EXPERIMENTS[args['experiment']]["script"](**args)
        

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run an experiment.')
    parser.add_argument('-e','--experiment', type=str, help='path to a python script')
    parser.add_argument('-d', '--data_path', type=str, help='path to a netcdf4 data')
    parser.add_argument('-m', '--metrics', nargs='+', help='names of metrics to use')
    parser.add_argument('-s', '--subset', help='should data be subset', action='store_true')
    parser.add_argument('-l', '--ls', help='list all experiments', action='store_true')
    args = parser.parse_args()
    args = vars(args)
    
    # list all experiements
    if args['ls']:
        list_all_experiments()

    # run experiment
    run_experiments(args)


