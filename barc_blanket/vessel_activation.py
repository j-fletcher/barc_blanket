import openmc
import openmc.stats
import openmc.deplete
import numpy as np

import os
import pickle

# Heavily based on John's stuff here: https://github.com/jlball/arc-nonproliferation/tree/master/openmc-scripts/arc-1/independent_depletion

def make_model():
    # Just quickfast snagging what Joe has to make a model

    if os.path.exists('model.xml'):
        model = openmc.Model.from_model_xml()
    else:
        ### MATERIALS ###

        # Vanadium alloy
        vanadium44 = openmc.Material(name='vanadium44')
        vanadium44.depletable = True
        vanadium44.add_element('V', 92.0, 'wo')
        vanadium44.add_element('Ti', 4.0, 'wo')
        vanadium44.add_element('Cr', 4.0, 'wo')
        vanadium44.set_density('g/cm3', 6.11)

        materials = openmc.Materials([vanadium44])

        major_radius = 680
        vv_thickness = 2

        # Set up our simplified reactor
        plasma_surface = openmc.Sphere(r=major_radius)
        vv_surface = openmc.Sphere(r=major_radius+vv_thickness)
        atmosphere_surface = openmc.Sphere(r=1000, boundary_type="vacuum")

        plasma_space = -plasma_surface
        vv_space = +plasma_surface & -vv_surface
        atmosphere_space = +vv_surface & -atmosphere_surface

        plasma_cell = openmc.Cell(name='plasma_cell', region=plasma_space, fill=None)
        vv_cell = openmc.Cell(name='vv_cell', region=vv_space, fill=vanadium44)
        atmosphere_cell = openmc.Cell(name='atmosphere_cell', region=atmosphere_space, fill=None)

        universe = openmc.Universe(cells=[plasma_cell, vv_cell, atmosphere_cell])
        geometry = openmc.Geometry(universe)

        # source definition
        point = openmc.stats.Point((0, 0, 0))
        source = openmc.IndependentSource(space=point)
        source.particle = 'neutron'
        source.angle = openmc.stats.Isotropic()
        source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

        # settings' settings
        settings = openmc.Settings(run_mode='fixed source')
        settings.photon_transport = False
        settings.source = source
        settings.batches = 40
        settings.particles = int(1e2) # modify this to shorten simulation, default was 1e6 

        model = openmc.Model(materials=materials, geometry=geometry,
                            settings=settings)

        model.export_to_model_xml()

        # Also save materials to file
        materials.export_to_xml()

    return model

def run_independent_vessel_activation(model:openmc.Model, days=365, num_timesteps=50, source_rate=2e20):
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
        The source rate of neutrons in the model. Default is 2e20 (for a 500MW fusion reactor)
    """

    openmc.config['cross_sections'] = '/home/zkeith/openmc_resources/endfb-viii.0-hdf5/cross_sections.xml'
    openmc.config['chain_file'] = '/home/zkeith/openmc_resources/chain_endfb80_sfr.xml'

    # Obtain a pointer to the vacuum vessel cell
    # TODO: there's definitely an easier way to do this
    vv_cell = next(iter(model._cells_by_name["vv_cell"]))

    # Check if flux and microscopic cross sections are present.
    # If not, calculate them
    # Otherwise, load them from file
    vv_flux_file = 'vv_flux.pkl'
    vv_microxs_file = 'vv_microxs.pkl'
    if os.path.exists(vv_flux_file) and os.path.exists(vv_microxs_file):
        with open(vv_flux_file, 'rb') as f:
            vv_flux = pickle.load(f)
        with open(vv_microxs_file, 'rb') as f:
            vv_microxs = pickle.load(f)
    else:
        vv_flux, vv_microxs = openmc.deplete.get_microxs_and_flux(model, [vv_cell])
        with open(vv_flux_file, 'wb') as f:
            pickle.dump(vv_flux, f)
        with open(vv_microxs_file, 'wb') as f:
            pickle.dump(vv_microxs, f)

    # TODO: programmatically calculate volume of cells
    # At the moment assuming R = 680 cm and a1 = 125 cm, a2 = 127 cm
    #vv_cell.fill.volume = (2*np.pi*680) * (np.pi*127**2 - np.pi*125**2)
    vv_cell.fill.volume = (4/3) * np.pi * (682**3 - 680**3)

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

def extract_activities():
    # Another thing taken from John: https://github.com/jlball/arc-nonproliferation/commit/04de395e19fd30344d9e5b2366918e149593b5d0
    openmc.config['cross_sections'] = '/home/zkeith/openmc_resources/endfb-viii.0-hdf5/cross_sections.xml'
    openmc.config['chain_file'] = '/home/zkeith/openmc_resources/chain_endfb80_sfr.xml'
    
    material = '0' # TODO: figure out how to get this properly

    # load results
    results = openmc.deplete.ResultsList.from_hdf5("depletion_results.h5")
    
    timesteps = results.get_times()
    activities = np.empty(len(timesteps))
    dists = []

    for i, step in enumerate(timesteps):
        materials = results.export_to_materials(i)

        vv_material = materials[0]

        activities[i] = vv_material.get_activity()
        dists.append(vv_material.get_decay_photon_energy())

    #activities = results.get_activity(material)

    return timesteps, activities, dists