from pathlib import Path
import openmc
import openmc.deplete
import numpy as np
from models.barc_model_simple_toroidal import make_model
from materials.blanket_depletion import gw_to_neutron_rate
from utilities import CROSS_SECTIONS, CHAIN_FILE

openmc.config['cross_sections'] = CROSS_SECTIONS
openmc.config['chain_file'] = CHAIN_FILE

# create model from barc_model_simple_toroidal.py

model = make_model(new_model_config={'particles':200,
                   'batches': 20})

# Set timesteps and source rates
# must have one source rate per timestep

timesteps_years = np.linspace(0, 200, 11)  # years

timesteps = np.array(timesteps_years) * 365 # convert to days

fusion_power = 2.2 # 2.2 GW
neutron_rate = gw_to_neutron_rate(fusion_power)
source_rates = np.ones_like(timesteps)*neutron_rate

# Setup CoupledOperator class 

op = openmc.deplete.CoupledOperator(model, 
                                    reduce_chain=True, 
                                    reduce_chain_level=3, 
                                    normalization_mode='source-rate')

# Set output directory 

output_dir = Path('./depletion_results')
op.output_dir=output_dir

# Run model

openmc.deplete.CECMIntegrator(op, timesteps, source_rates=source_rates, timestep_units='d').integrate()