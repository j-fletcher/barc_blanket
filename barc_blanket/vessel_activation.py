import openmc
import openmc.stats
import openmc.deplete
import numpy as np

import os

from barc_blanket.utilities import CROSS_SECTIONS, CHAIN_FILE

# Heavily based on John's stuff here: https://github.com/jlball/arc-nonproliferation/tree/master/openmc-scripts/arc-1/independent_depletion

def run_independent_vessel_activation(model:openmc.Model, days=365, num_timesteps=50, source_rate=3.6e20):
    """ Run the vessel activation after a certain number of days.

    Parameters:
    -----------
    model : openmc.Model
        The model to get the vessel activation from.
    days : int
        The number of days to run the model for.
    num_timesteps : int
        The number of timesteps to run the model for.
    source_rate : float
        The source rate of neutrons in the model. Default is 3.6e20 (for 1 GW fusion power)
    """

    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    # Obtain a pointer to the vacuum vessel cell
    # TODO: there's definitely an easier way to do this
    vv_cell = next(iter(model._cells_by_name["inboard_vv_cell"]))

    # Check if flux and microscopic cross sections are present.
    # If not, calculate them
    # Otherwise, load them from file
    vv_flux_file = 'vv_cell.npy'
    vv_microxs_file = 'vv_cell_microxs.csv'
    if os.path.exists(vv_flux_file) and os.path.exists(vv_microxs_file):
        with open(vv_flux_file, 'rb') as f:
            # pickle.load(f)
            vv_flux = [np.load(vv_flux_file)]#pickle.load(f)
        with open(vv_microxs_file, 'rb') as f:
            #vv_microxs = pickle.load(f)
            vv_microxs = [openmc.deplete.MicroXS.from_csv(vv_microxs_file)]
    else:
        vv_flux, vv_microxs = openmc.deplete.get_microxs_and_flux(model, [vv_cell])
        np.save(vv_flux_file, vv_flux[0])
        vv_microxs[0].to_csv(vv_microxs_file)

    # TODO: programmatically calculate volume of cells
    # At the moment assuming R = 680 cm and a1 = 125 cm, a2 = 127 cm
    #vv_cell.fill.volume = (2*np.pi*680) * (np.pi*127**2 - np.pi*125**2)
    vv_cell.fill.volume = (4/3) * np.pi * (682**3 - 680**3)
    print(vv_cell.fill)

    # Perform depletion (CHECK NORMALIZATION MODE)
    vv_operator = openmc.deplete.IndependentOperator(openmc.Materials([vv_cell.fill]),
                                                    vv_flux,
                                                    vv_microxs,
                                                    normalization_mode='source-rate',
                                                    reduce_chain=True,
                                                    reduce_chain_level=5) # TODO: figure out what this does and why we set to 5
    
    time_steps = [days/num_timesteps] * num_timesteps
    source_rates = np.ones(num_timesteps) * source_rate

    vv_integrator = openmc.deplete.PredictorIntegrator(vv_operator, 
                                                       time_steps,
                                                       source_rates=source_rates,
                                                       timestep_units='d')
    
    vv_integrator.integrate()

def extract_activities(model:openmc.Model):
    # Get the total activity in the vacuum vessel cell from depletion results
    # Another thing taken from John: https://github.com/jlball/arc-nonproliferation/commit/04de395e19fd30344d9e5b2366918e149593b5d0
    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    # load results
    results = openmc.deplete.Results("depletion_results.h5")
    
    timesteps = results.get_times()
    activities = np.empty(len(timesteps))
    

    vv_cell = next(iter(model._cells_by_name["vv_cell"]))
    print("second fill")
    print(vv_cell.fill)
    activities = results.get_activity(vv_cell.fill)

    return timesteps, activities

def extract_decay_photon_energies():
    # Get the decay photon energies from the depletion results
    
    results = openmc.deplete.Results("depletion_results.h5")
    timesteps = results.get_times()
    activities = np.empty(len(timesteps))

    dists = []

    for i, step in enumerate(timesteps):
        materials = results.export_to_materials(i)

        vv_material = materials[0]

        activities[i] = vv_material.get_activity()
        dists.append(vv_material.get_decay_photon_energy())

    return timesteps, dists

def extract_nuclides(model:openmc.Model):
    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    results = openmc.deplete.Results("depletion_results.h5")
    timesteps = results.get_times()
    vv_cell = next(iter(model._cells_by_name["vv_cell"]))
    
    nuc_atoms = {}
    for nuclide_tuple in vv_cell.fill.nuclides:
        nuclide_name = nuclide_tuple[0]
        nuc_atoms[nuclide_name] = results.get_atoms(vv_cell.fill, nuc=nuclide_name)[1]

    return timesteps, nuc_atoms
