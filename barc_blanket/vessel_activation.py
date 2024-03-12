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
    model_config : dict
        The configuration of the model, detailing dimensions of geometry
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
    vv_cell = next(iter(model._cells_by_name["vv_cell"]))
    bv_cell = next(iter(model._cells_by_name["bv_cell"]))

    # Check if flux and microscopic cross sections are present.
    # If not, calculate them
    # Otherwise, load them from file
    fluxes_file = 'fluxes.npy'
    # TODO: should be able to programmatically put all the microxs in one file, but for now we'll just do it separately
    vv_microxs_file = 'vv_microxs.csv'
    bv_microxs_file = 'bv_microxs.csv'
    if os.path.exists(fluxes_file) and os.path.exists(vv_microxs_file) and os.path.exists(bv_microxs_file):
        with open(fluxes_file, 'rb') as f:
            fluxes = np.load(fluxes_file)
        with open(vv_microxs_file, 'rb') as f:
            vv_microxs = openmc.deplete.MicroXS.from_csv(vv_microxs_file)
        with open(bv_microxs_file, 'rb') as f:
            bv_microxs = openmc.deplete.MicroXS.from_csv(bv_microxs_file)
    else:
        fluxes, microxs = openmc.deplete.get_microxs_and_flux(model, [vv_cell, bv_cell])
        np.save(fluxes_file, fluxes)
        vv_microxs = microxs[0]
        bv_microxs = microxs[1]
        vv_microxs.to_csv(vv_microxs_file)
        bv_microxs.to_csv(bv_microxs_file)

    # Perform depletion (CHECK NORMALIZATION MODE)
    operator = openmc.deplete.IndependentOperator(openmc.Materials([vv_cell.fill, bv_cell.fill]),
                                                    fluxes,
                                                    [vv_microxs, bv_microxs],
                                                    normalization_mode='source-rate',
                                                    reduce_chain=True,
                                                    reduce_chain_level=5) # TODO: figure out what this does and why we set to 5
    
    time_steps = [days/num_timesteps] * num_timesteps
    source_rates = np.ones(num_timesteps) * source_rate

    integrator = openmc.deplete.PredictorIntegrator(operator, 
                                                       time_steps,
                                                       source_rates=source_rates,
                                                       timestep_units='d')
    
    integrator.integrate()

def extract_activities(model:openmc.Model, cell_name:str="vv_cell"):
    # Get the total activity from a specified cell
    # Another thing taken from John: https://github.com/jlball/arc-nonproliferation/commit/04de395e19fd30344d9e5b2366918e149593b5d0
    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    # load results
    results = openmc.deplete.Results("depletion_results.h5")
    
    timesteps = results.get_times()
    activities = np.empty(len(timesteps))
    

    cell = next(iter(model._cells_by_name[cell_name]))
    activities = results.get_activity(cell.fill)

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


def extract_activated_model():
    # Make a new model from the activated materials



def run_independent_vessel_decay(model:openmc.Model, days=365, num_timesteps=50):
    # Run the decay of the vessel after a certain number of days

    