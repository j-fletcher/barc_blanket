from pathlib import Path
import openmc
import openmc.deplete
import numpy as np
from barc_blanket.models.barc_model_final import make_model
from barc_blanket.utilities import CROSS_SECTIONS, CHAIN_FILE
from barc_blanket.models.materials_tools import burner_mixture

openmc.config['cross_sections'] = CROSS_SECTIONS
openmc.config['chain_file'] = CHAIN_FILE

# create model from barc_model_simple_toroidal.py

model = make_model(new_model_config={
    'particles': int(1e4),
    'blanket_material': burner_mixture(0.1, id=6),
    'removed_U238': 1.0
    })

# Set timesteps and source rates
# must have one source rate per timestep
nt = 20
t_max = 100 # years
timesteps_years = t_max * np.ones(nt) / nt  # years
timesteps = np.array(timesteps_years) * 365 # convert to days
print(timesteps)
print(sum(timesteps_years))

efus = 17.6e6  # eV
ev2j = 1.60218e-19
Pfus = 2200e6  # W
neutron_rate = Pfus / (efus * ev2j)  # n/s

source_rates = np.ones_like(timesteps)*neutron_rate

# Setup CoupledOperator class 

op = openmc.deplete.CoupledOperator(model, 
                                    reduce_chain=True, 
                                    reduce_chain_level=3, 
                                    normalization_mode='source-rate')

# Set output directory 

output_dir = Path('./depletion_results_final/barc-v7-10-percent-waste-sludge-U238-removed')
op.output_dir=output_dir

# Run model

integrator = openmc.deplete.CECMIntegrator(op, timesteps, source_rates=source_rates, timestep_units='d')

integrator.integrate()