import os
from contextlib import contextmanager

CROSS_SECTIONS = '/usr/local/share/xs_data/endfb80_hdf5/cross_sections.xml'
CHAIN_FILE = '/home/zkeith/openmc_resources/chain_endfb80_sfr.xml'

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)