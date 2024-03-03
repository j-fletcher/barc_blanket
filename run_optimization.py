import os
import sys
import yaml
import optuna

from barc_blanket.optimize_geometry import make_model, evaluate_metric

def objective(trial, sweep_config):
    """ Objective function for the optimization

    Parameters:
    ----------
    trial : optuna.Trial
        An optuna trial object
    sweep_config : dict
        A dictionary containing the sweep configuration

    Returns:
    -------
    metric_val: float
        The value of the metric calculated for the parameters in this trial
    """

    # Obtain the values of parameters from the trial
    model_config = {}
    parameters = sweep_config['parameters']
    parameter_names = list(parameters.keys())

    for parameter_name in parameter_names:
        parameter = parameters[parameter_name]
        distribution_type = parameter['distribution']

        if distribution_type == "int":
            min = parameter['min']
            max = parameter['max']
            chosen_value = trial.suggest_int(parameter_name, min, max)
        elif distribution_type == "float":
            min = parameter['min']
            max = parameter['max']
            log = parameter['log']
            chosen_value = trial.suggest_float(parameter_name, min, max, log=log)
        elif distribution_type == "categorical":
            values = parameter['values']
            chosen_value = trial.suggest_categorical(parameter_name, values)
        else:
            raise ValueError(f"Invalid distribution type: {distribution_type}")
        
        model_config[parameter_name] = chosen_value

    # Create the model and evaluate the metric
    try:
        model = make_model(model_config)
        metric_val = evaluate_metric(model, sweep_config['metric'])
    except MemoryError as e:
        print(f"Ran out of memory for trial {trial.number}")
        print(e)
        metric_val = float('nan')
    except Exception as e:
        print(f"Error in trial {trial.number}, pruning...")
        print(e)
        # If anything oges wrong during training or validation, say that the trial was pruned
        # This should make Optuna try a different set of parameters to avoid errors
        raise optuna.TrialPruned()
    
    return metric_val

def main(sweep_config_path, num_trials=1):
    """Run a parameter sweep using Optuna

    Parameters:
    ----------
    sweep_config_path : str
        Path to the sweep configuration file
        This configuration file should be a YAML file with the following structure:
        ```
        metric: metric name
        direction: minimize or maximize
        parameters:
          parameter1:
            distribution: int
            max: 100
            min: 0
          parameter2:
            distribution: float
            log: true for log scale, false for linear scale
            max: 1.0
            min: 0.0
          parameter3:
            distribution: categorical
            values:
            - value1
            - value2
        ```
    num_trials:int, default=1
        Number of trials to run. This will add num_trials to the existing trials in the sweep_results.db
    """

    # Load the config
    sweep_config = yaml.safe_load(open(sweep_config_path, "r"))

    # Create storage for trial results that can support concurrent writes
    sweep_results_path = f"{sweep_config_path}/sweep_results.db"
    lock_obj = optuna.storages.JournalFileOpenLock(sweep_results_path)
    storage = optuna.storages.JournalStorage(
        optuna.storages.JournalFileStorage(sweep_results_path, lock_obj=lock_obj)
    )

    study = optuna.create_study(
        storage=storage, 
        study_name=f"{sweep_config_path}",
        direction=sweep_config['direction'],
        load_if_exists=True
    )

    study.optimize(lambda trial: objective(trial, sweep_config), n_trials=num_trials)


if __name__ == "__main__":
    if len(sys.argv) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], int(sys.argv[2]))
    else:
        print("Usage: python run_optimization.py <sweep_config_path> <num_trials>")
        sys.exit(1)