import os
from contextlib import contextmanager

CROSS_SECTIONS = '/home/hallj/endfb80_hdf5/cross_sections.xml'
CHAIN_FILE = '/home/hallj/barc_blanket/TENDL_cross_sections/chain_endfb80_sfr_barc.xml'

@contextmanager
def working_directory(directory):
    owd = os.getcwd()
    try:
        os.chdir(directory)
        yield directory
    finally:
        os.chdir(owd)