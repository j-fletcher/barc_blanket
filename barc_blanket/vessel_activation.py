import openmc
import openmc.stats
import openmc.deplete
import numpy as np
import tempfile
from typing import Optional, List
from barc_blanket.utilities import change_directory

import os

from barc_blanket.utilities import CROSS_SECTIONS, CHAIN_FILE

from openmc_regular_mesh_plotter import plot_mesh_tally
from matplotlib.colors import LogNorm

# Heavily based on John's stuff here: https://github.com/jlball/arc-nonproliferation/tree/master/openmc-scripts/arc-1/independent_depletion

def run_independent_vessel_activation(model:openmc.Model, days=365, num_timesteps=50, times=None, source_rate=3.6e20):
    """ Run the vessel activation after a certain number of days.

    Parameters:
    -----------
    model : openmc.Model
        The model to get the vessel activation from.
    days : int
        The number of days to run the model for. Default is 365.
    num_timesteps : int
        The number of timesteps to run the model for. Default is 50.
    times : list
        The times to evaluate activation. Default is None. If not none, it will override days and num_timesteps.
    source_rate : float
        The source rate of neutrons in the model. Default is 3.6e20 (for 1 GW fusion power)
    """

    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    # Obtain a pointer to the vacuum vessel cell
    first_wall_cell = next(iter(model._cells_by_name["first_wall_cell"]))
    vacuum_vessel_cell = next(iter(model._cells_by_name["vacuum_vessel_cell"]))
    blanket_vessel_cell = next(iter(model._cells_by_name["blanket_vessel_cell"]))

    # Check if flux and microscopic cross sections are present.
    # If not, calculate them
    # Otherwise, load them from file
    fluxes_file = 'fluxes.npy'
    # TODO: should be able to programmatically put all the microxs in one file, but for now we'll just do it separately
    first_wall_file = 'first_wall_microxs.csv'
    vacuum_vessel_microxs_file = 'vacuum_vessel_microxs.csv'
    blanket_vessel_microxs_file = 'blanket_vessel_microxs.csv'
    if os.path.exists(fluxes_file) and os.path.exists(first_wall_file) and os.path.exists(vacuum_vessel_microxs_file) and os.path.exists(blanket_vessel_microxs_file):
        with open(fluxes_file, 'rb') as f:
            fluxes = np.load(fluxes_file)
        with open(first_wall_file, 'rb') as f:
            first_wall_microxs = openmc.deplete.MicroXS.from_csv(first_wall_file)
        with open(vacuum_vessel_microxs_file, 'rb') as f:
            vacuum_vessel_microxs = openmc.deplete.MicroXS.from_csv(vacuum_vessel_microxs_file)
        with open(blanket_vessel_microxs_file, 'rb') as f:
            blanket_vessel_microxs = openmc.deplete.MicroXS.from_csv(blanket_vessel_microxs_file)
    else:
        fluxes, microxs = openmc.deplete.get_microxs_and_flux(model, [first_wall_cell, vacuum_vessel_cell, blanket_vessel_cell])
        np.save(fluxes_file, fluxes)
        first_wall_microxs = microxs[0]
        vacuum_vessel_microxs = microxs[1]
        blanket_vessel_microxs = microxs[2]
        first_wall_microxs.to_csv(first_wall_file)
        vacuum_vessel_microxs.to_csv(vacuum_vessel_microxs_file)
        blanket_vessel_microxs.to_csv(blanket_vessel_microxs_file)

    # Perform depletion (CHECK NORMALIZATION MODE)
    operator = openmc.deplete.IndependentOperator(openmc.Materials([first_wall_cell.fill, vacuum_vessel_cell.fill, blanket_vessel_cell.fill]),
                                                    fluxes,
                                                    [first_wall_microxs, vacuum_vessel_microxs, blanket_vessel_microxs],
                                                    normalization_mode='source-rate',
                                                    reduce_chain=True,
                                                    reduce_chain_level=5) # TODO: figure out what this does and why we set to 5
    
    # 'timestep' is the actual time of the depletion step
    # 'timediff' is the difference in time between the current and previous depletion step
    if times is None:
        timesteps = [days/num_timesteps] * num_timesteps
    else:
        timesteps = np.diff(times)
    source_rates = np.ones(len(timesteps)) * source_rate
    
    integrator = openmc.deplete.PredictorIntegrator(operator, 
                                                       timesteps,
                                                       source_rates=source_rates,
                                                       timestep_units='d')
    
    integrator.integrate()

def run_independent_vessel_decay(model:openmc.Model, results, days=365, num_timesteps=50, times=None):
    """ Run the vessel decay after a certain number of days.

    Parameters:
    -----------
    model : openmc.Model
        The model to get the vessel activation from.
    days : int
        The number of days to run the model for.
    num_timesteps : int
        The number of timesteps to run the model for.
    times : list
        The times to evaluate activation. If not none, it will override days and num_timesteps.
    """

    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    # Obtain a pointer to the vessel cells
    first_wall_cell = next(iter(model._cells_by_name["first_wall_cell"]))
    vacuum_vessel_cell = next(iter(model._cells_by_name["vacuum_vessel_cell"]))
    blanket_vessel_cell = next(iter(model._cells_by_name["blanket_vessel_cell"]))

    # This is a hack until I can figure out how to make it actually work
    fluxes_file = '../independent_vessel_activation/fluxes.npy'
    # TODO: should be able to programmatically put all the microxs in one file, but for now we'll just do it separately
    first_wall_file = '../independent_vessel_activation/first_wall_microxs.csv'
    vacuum_vessel_microxs_file = '../independent_vessel_activation/vacuum_vessel_microxs.csv'
    blanket_vessel_microxs_file = '../independent_vessel_activation/blanket_vessel_microxs.csv'
    if os.path.exists(fluxes_file) and os.path.exists(first_wall_file) and os.path.exists(vacuum_vessel_microxs_file) and os.path.exists(blanket_vessel_microxs_file):
        with open(fluxes_file, 'rb') as f:
            fluxes = np.load(fluxes_file)
        with open(first_wall_file, 'rb') as f:
            first_wall_microxs = openmc.deplete.MicroXS.from_csv(first_wall_file)
        with open(vacuum_vessel_microxs_file, 'rb') as f:
            vacuum_vessel_microxs = openmc.deplete.MicroXS.from_csv(vacuum_vessel_microxs_file)
        with open(blanket_vessel_microxs_file, 'rb') as f:
            blanket_vessel_microxs = openmc.deplete.MicroXS.from_csv(blanket_vessel_microxs_file)

    fluxes = np.zeros(fluxes.shape)

    # Perform depletion (CHECK NORMALIZATION MODE)
    operator = openmc.deplete.IndependentOperator(openmc.Materials([first_wall_cell.fill, vacuum_vessel_cell.fill, blanket_vessel_cell.fill]),
                                                    fluxes,
                                                    [first_wall_microxs, vacuum_vessel_microxs, blanket_vessel_microxs],
                                                    normalization_mode='source-rate',
                                                    prev_results=results,
                                                    reduce_chain=True,
                                                    reduce_chain_level=5) # TODO: figure out what this does and why we set to 5
    
    if times is None:
        timesteps = [days/num_timesteps] * num_timesteps
    else:
        timesteps = np.diff(times)
    source_rates = np.ones(len(timesteps)) # source rate is very close to 0 so decay can happen

    integrator = openmc.deplete.PredictorIntegrator(operator, 
                                                       timesteps,
                                                       source_rates=source_rates,
                                                       timestep_units='d')
    
    integrator.integrate()

def extract_activities(model:openmc.Model, cell_name:str="blanket_vessel_cell"):
    # Get the total activity from a specified cell
    # Another thing taken from John: https://github.com/jlball/arc-nonproliferation/commit/04de395e19fd30344d9e5b2366918e149593b5d0
    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    # load results
    results = openmc.deplete.Results("depletion_results.h5")
    
    times = results.get_times()
    activities = np.empty(len(times))

    cell = next(iter(model._cells_by_name[cell_name]))
    activities = results.get_activity(cell.fill)

    return times, activities

def extract_decay_heat(model:openmc.Model, cell_name:str="blanket_vessel_cell"):
    """ Get the decay heat from a specified cell.
    
    Parameters:
    -----------
    model : openmc.Model
        The model to get the decay heat from.
    cell_name : str
        The name of the cell to get the decay heat from.

    Returns:
    --------
    times : list
        The times at which the decay heat was calculated.
        Starts at 0, which is either the initial activation time or the initial decay time.
    decay_heats : list
        The decay heat at each time, in Watts.
    """

    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE
    
    results = openmc.deplete.Results("depletion_results.h5")

    cell = next(iter(model._cells_by_name[cell_name]))
    decay_heat_array = results.get_decay_heat(cell.fill)

    # Limit real times to only where decay_heat_array[1] is not nan
    times = []
    decay_heats = []
    for i, heat in enumerate(decay_heat_array[1]):
        if not np.isnan(decay_heat_array[1][i]):
            times.append(decay_heat_array[0][i])
            decay_heats.append(decay_heat_array[1][i])
    
    times = times - times[0]  # start at 0
    return times, decay_heats

def extract_decay_photon_energies():
    # Get the decay photon energies from the depletion results
    
    results = openmc.deplete.Results("depletion_results.h5")
    times = results.get_times()
    activities = np.empty(len(times))

    dists = []

    for i, time in enumerate(times):
        materials = results.export_to_materials(i)

        vacuum_vessel_material = materials[0]

        activities[i] = vacuum_vessel_material.get_activity()
        dists.append(vacuum_vessel_material.get_decay_photon_energy())

    return times, dists

def extract_original_nuclides(model:openmc.Model, cell_name:str="blanket_vessel_cell"):
    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    results = openmc.deplete.Results("depletion_results.h5")
    times = results.get_times()
    cell = next(iter(model._cells_by_name[cell_name]))
    
    nuc_atoms = {}
    for nuclide_tuple in cell.fill.nuclides:
        nuclide_name = nuclide_tuple[0]
        nuc_atoms[nuclide_name] = results.get_atoms(cell.fill, nuc=nuclide_name)[1]

    return times, nuc_atoms

def extract_nuclides(model:openmc.Model, cell_name:str="blanket_vessel_cell", nuclide_names:list=["H-3", "He-4"]):
    openmc.config['cross_sections'] = CROSS_SECTIONS
    openmc.config['chain_file'] = CHAIN_FILE

    results = openmc.deplete.Results("depletion_results.h5")
    times = results.get_times()
    cell = next(iter(model._cells_by_name[cell_name]))
    
    nuc_atoms = {}
    for nuclide_name in nuclide_names:
        nuc_atoms[nuclide_name] = np.empty(len(times))
        nuc_atoms[nuclide_name] = results.get_atoms(cell.fill, nuc=nuclide_name)[1]

    return times, nuc_atoms

def plot_2d_dose(statepoint, mesh):
    photon_tally = statepoint.get_tally(name="photon_dose_on_mesh")

    # normalising this tally is a little different to other examples as the source strength has been using units of photons per second.
    # tally.mean is in units of pSv-cm3/source photon.
    # as source strength is in photons_per_second this changes units to pSv-/second

    # multiplication by pico_to_micro converts from (pico) pSv/s to (micro) uSv/s
    # dividing by mesh voxel volume cancles out the cm3 units
    # could do the mesh volume scaling on the plot and vtk functions but doing it here instead
    pico_to_micro = 1e-6
    seconds_to_hours = 60*60
    scaling_factor = (seconds_to_hours * pico_to_micro) / mesh.volumes[0][0][0]

    plot = plot_mesh_tally(
            tally=photon_tally,
            basis="xz",
            # score='flux', # only one tally so can make use of default here
            value="mean",
            colorbar_kwargs={
                'label': "Decay photon dose [µSv/h]",
            },
            norm=LogNorm(),
            volume_normalization=False,  # this is done in the scaling_factor
            scaling_factor=scaling_factor,
        )
    plot.figure.savefig(f'shut_down_dose_map_timestep_{i_cool}')

def get_cell_volume_in_mesh(
            self,
            model: openmc.Model,
            n_samples: int = 10_000,
            prn_seed: Optional[int] = None,
            **kwargs
    ) -> List[openmc.Material]:
        """ Get the volume of a particular cell 
        Based on the code here:
        https://github.com/openmc-dev/openmc/pull/2971/files#diff-967783d59b58404de3b672e391edb35b65e00bf29f89f11edfa8b44469494c50
        
        Generate homogenized materials over each element in a mesh.
        .. versionadded:: 0.14.1
        Parameters
        ----------
        model : openmc.Model
            Model containing materials to be homogenized and the associated
            geometry.
        n_samples : int
            Number of samples in each mesh element.
        prn_seed : int, optional
            Pseudorandom number generator (PRNG) seed; if None, one will be
            generated randomly.
        **kwargs
            Keyword-arguments passed to :func:`openmc.lib.init`.
        Returns
        -------
        material_volumes: dict
            A dictionary where the keys are the mesh element IDs and the values
            are the homogenized materials over the element.


        """
        import openmc.lib

        with change_directory(tmpdir=True):
            # In order to get mesh into model, we temporarily replace the
            # tallies with a single mesh tally using the current mesh
            original_tallies = model.tallies
            new_tally = openmc.Tally()
            new_tally.filters = [openmc.MeshFilter(self)]
            new_tally.scores = ['flux']
            model.tallies = [new_tally]

            # Export model to XML
            model.export_to_model_xml()

            # Get material volume fractions
            openmc.lib.init(**kwargs)
            mesh = openmc.lib.tallies[new_tally.id].filters[0].mesh
            mat_volume_by_element = [
                [
                    (mat.id if mat is not None else None, volume)
                    for mat, volume in mat_volume_list
                ]
                for mat_volume_list in mesh.material_volumes(n_samples, prn_seed)
            ]
            openmc.lib.finalize()

            # Restore original tallies
            model.tallies = original_tallies

        # Create homogenized material for each element
        materials = model.geometry.get_all_materials()
        homogenized_materials = []
        for mat_volume_list in mat_volume_by_element:
            material_ids, volumes = [list(x) for x in zip(*mat_volume_list)]
            total_volume = sum(volumes)

            # Check for void material and remove
            try:
                index_void = material_ids.index(None)
            except ValueError:
                pass
            else:
                material_ids.pop(index_void)
                volumes.pop(index_void)

            # Compute volume fractions
            volume_fracs = np.array(volumes) / total_volume

            # Get list of materials and mix 'em up!
            mats = [materials[uid] for uid in material_ids]
            homogenized_mat = openmc.Material.mix_materials(
                mats, volume_fracs, 'vo'
            )
            homogenized_mat.volume = total_volume
            homogenized_materials.append(homogenized_mat)

        return homogenized_materials