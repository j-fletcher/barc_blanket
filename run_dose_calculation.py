import os
import openmc
import matplotlib.pyplot as plt

from barc_blanket.vessel_activation import extract_decay_photon_energies

# Create a place to put all the files we'll be working with for depletion
# TODO: Implement with context switcher
working_directory = "depletion"
if not os.path.exists(working_directory):
    os.makedirs(working_directory)
os.chdir(working_directory)

# Load model
model = openmc.model.Model.from_xml("model.xml")

timesteps, dists = extract_decay_photon_energies(model)





