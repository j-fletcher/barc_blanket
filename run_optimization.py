import sys
import yaml
import optuna

def main(sweep_config_path):


    # Load the sweep configuration file
    sweep_config = yaml.safe_load(open(sweep_config_path, 'r'))


    lock_obj = optuna.storages.JournalFileOpenLock(database_path)


if __name__ == "__main__":
    main(sys.argv[1])