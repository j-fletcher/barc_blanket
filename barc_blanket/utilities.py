import os
from contextlib import contextmanager

CROSS_SECTIONS = '/home/hallj/endfb71_hdf5/cross_sections.xml'
CHAIN_FILE = '/usr/local/share/xs_data/depletion_chains/chain_endfb80_pwr.xml'

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)