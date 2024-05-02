from pathlib import Path
import openmc
import openmc.deplete
import numpy as np
from barc_blanket.models.barc_model_simple_toroidal import make_model
from barc_blanket.utilities import CROSS_SECTIONS, CHAIN_FILE

openmc.config['cross_sections'] = CROSS_SECTIONS
openmc.config['chain_file'] = CHAIN_FILE

# create model from barc_model_simple_toroidal.py

model = make_model(new_model_config={
    'slurry_ratio': 0.05,
    'removed_U': 0.0,
    'removed_Pu': 0.0,
    'particles': int(20)
    })

blanket_cell = model.geometry.get_cells_by_name("blanket_cell")[0]
print(blanket_cell.fill.volume)
cooling_channel_cell = model.geometry.get_cells_by_name("cooling_channel_cell")[0]

flux_list, microxs_list = openmc.deplete.get_microxs_and_flux(
    model, 
    domains=[blanket_cell],
    chain_file='/home/hallj/barc_blanket/TENDL_cross_sections/chain_endfb80_sfr_barc.xml'
    )

flux = flux_list[0] # n-cm / src
efus = 17.6e6  # eV
ev2j = 1.60218e-19
Pfus = 2200e6  # W
neutron_rate = Pfus / (efus / ev2j)  # n/s

scaled_flux = flux / blanket_cell.fill.volume * neutron_rate
print(scaled_flux)

for i in microxs_list:
    i.to_csv("microxs.csv")