import openmc.deplete
from barc_blanket.utilities import working_directory

with working_directory("independent_vessel_decay"):
    # 1. Get the activated materials from the depletion results
    
    activation_results = openmc.deplete.Results("../independent_vessel_activation/depletion_results.h5")

    # Get model from the last timestep of depletion results
    

#2. Make a new model with the activated materials

#3. Run the new model with the activated materials
