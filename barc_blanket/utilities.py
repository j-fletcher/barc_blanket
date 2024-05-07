import os
from contextlib import contextmanager

CROSS_SECTIONS = os.environ['OPENMC_CROSS_SECTIONS']
CHAIN_FILE = os.environ['OPENMC_CHAIN_FILE']

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)