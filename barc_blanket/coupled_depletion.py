from pathlib import Path
import openmc
import openmc.deplete
import numpy as np
from models.barc_model_simple_toroidal import make_model
from utilities import CROSS_SECTIONS, CHAIN_FILE

openmc.config['cross_sections'] = CROSS_SECTIONS
openmc.config['chain_file'] = CHAIN_FILE

# create model from barc_model_simple_toroidal.py

model = make_model()

# Set timesteps and source rates
# must have one source rate per timestep

timesteps = [10.0, 10.0, 10.0]  # days

efus = 17.6e6  # eV
ev2j = 1.60218e-19
Pfus = 1000e6  # W
neutron_rate = Pfus / efus / ev2j  # n/s

source_rates = np.ones_like(timesteps)*neutron_rate

# Setup CoupledOperator class 

op = openmc.deplete.CoupledOperator(model, 
                                    reduce_chain=True, 
                                    reduce_chain_level=3, 
                                    normalization_mode='source-rate')

# Set output directory 

output_dir = Path('./depletion_results_tank')
op.output_dir=output_dir

# Run model

openmc.deplete.CECMIntegrator(op, timesteps, source_rates=source_rates, timestep_units='d').integrate()