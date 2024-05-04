import os
from contextlib import contextmanager

CROSS_SECTIONS = '/usr/local/share/xs_data/endfb80_hdf5/cross_sections.xml'
CHAIN_FILE = '/home/zkeith/proj/barc_blanket2/barc_blanket/TENDL_cross_sections/chain_endfb80_sfr_barc.xml'

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)