from pathlib import Path
import openmc
import openmc.deplete
import numpy as np
from barc_blanket.models.barc_model_simple_toroidal import make_model

# %%
loaded_model = make_model() #openmc.Model.from_model_xml('../models/model.xml')
# %%

timesteps = [10.0, 10.0, 10.0]  # days
#timesteps = np.logspace(0,7.49,5)

efus = 17.6e6  # eV
ev2j = 1.60218e-19
Pfus = 1000e6  # W

neutron_rate = Pfus / efus / ev2j  # n/s

source_rates = np.ones_like(timesteps)*neutron_rate

output_dir = Path('./depletion_results')

op = openmc.deplete.CoupledOperator(loaded_model, 
                                    reduce_chain=True, 
                                    reduce_chain_level=3, 
                                    normalization_mode='source-rate',
                                    output_dir=output_dir)
#openmc.config['chain_file'] = '/Users/hallj/barc_blanket/chain_endfb71_pwr.xml'
#os.environ["OMP_NUM_THREADS"] = '8'

openmc.deplete.CECMIntegrator(op, timesteps, source_rates=source_rates, timestep_units='d').integrate()
