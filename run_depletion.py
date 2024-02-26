from barc_blanket.vessel_activation import make_model
import os


# Create a place to put all the files we'll be working with for depletion (TODO: standardize this structure with the rest of the group)
working_directory = "depletion"
if not os.path.exists(working_directory):
    os.makedirs(working_directory)
os.chdir(working_directory)

make_model()