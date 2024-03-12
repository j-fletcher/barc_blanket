"""
MIT License

Copyright (c) 2021 fusion-energy neutronics-workshop contributors

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

https://github.com/fusion-energy/neutronics-workshop/blob/reshuffle/tasks/task_11_CSG_shut_down_dose_tallies/1_cell_based_shut_down_dose_rate_example.py
"""

import os
import openmc
import openmc.model
import openmc.deplete
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from barc_blanket.vessel_activation import CHAIN_FILE, CROSS_SECTIONS

openmc.config['cross_sections'] = CROSS_SECTIONS
openmc.config['chain_file'] = CHAIN_FILE

# Create a place to put all the files we'll be working with for depletion
# TODO: Implement with context switcher
working_directory = "depletion"
if not os.path.exists(working_directory):
    os.makedirs(working_directory)
os.chdir(working_directory)

# Load model
model = openmc.model.Model.from_model_xml("model.xml")

# Create decay gamma simulation
gamma_settings = openmc.Settings()
gamma_settings.particles = 1000
gamma_settings.batches = 100
gamma_settings.run_mode = "fixed source"

mesh = openmc.RegularMesh().from_domain(
    model.geometry,
    dimension=[20, 20, 10], # 100 voxels in each axis direction
)

# Add dose tally to the regular mesh
# AP, PA, LLAT, RLAT, ROT, ISO are ICRP incident dose field directions, AP is front facing
energies, pSv_cm2 = openmc.data.dose_coefficients(particle="photon", geometry="AP")
dose_filter = openmc.EnergyFunctionFilter(
    energies, pSv_cm2, interpolation="cubic" # Apparently cubic interpolation is recommended by ICRP
)
particle_filter = openmc.ParticleFilter(["photon"])
mesh_filter = openmc.MeshFilter(mesh)
flux_tally = openmc.Tally()
flux_tally.filters = [mesh_filter, dose_filter, particle_filter]#, particle_filter, mesh_filter]
flux_tally.scores = ["flux"]
flux_tally.name = "photon_dose_on_mesh"

tallies = openmc.Tallies([flux_tally])

activated_cell_ids = [c.id for c in model.geometry.get_all_material_cells().values() if c.fill.depletable]
cells = model.geometry.get_all_cells()
activated_cells = [cells[uid] for uid in activated_cell_ids]

results = openmc.deplete.Results("depletion_results.h5")
timesteps = results.get_times()

for i_cool in range(len(timesteps)-1, len(timesteps)):

    # range starts at 1 to skip the first step as that is an irradiation step and there is no
    # decay gamma source from the stable material at that time
    # also there are no decay products in this first timestep for this model

    photon_sources_for_timestep = []
    print(f"making photon source for timestep {i_cool}")

    all_activated_materials_in_timestep = []

    for activated_cell_id in activated_cell_ids:
        # gets the material id of the material filling the cell
        material_id = cells[activated_cell_id].fill.id

        # gets the activated material using the material id
        activated_mat = results[i_cool].get_material(str(material_id))
        # gets the energy and probabilities for the 
        energy = activated_mat.get_decay_photon_energy()
        strength = energy.integral()

        if strength > 0.:  # only makes sources for 
            space = openmc.stats.Box(*cells[activated_cell_id].bounding_box)
            source = openmc.IndependentSource(
                space=space,
                energy=energy,
                particle="photon",
                strength=strength,
                domains=[cells[activated_cell_id]],
            )
            photon_sources_for_timestep.append(source)

    gamma_settings.source = photon_sources_for_timestep

    # TODO: run with depleted material, not pristine material
    model_gamma = openmc.Model(model.geometry, model.materials, gamma_settings, tallies)

    model_gamma.run()

pico_to_micro = 1e-6
seconds_to_hours = 60*60

# You may wish to plot the dose tally on a mesh, this package makes it easy to include the geometry with the mesh tally
from openmc_regular_mesh_plotter import plot_mesh_tally
#for i_cool in range(1, len(timesteps)):
with openmc.StatePoint('statepoint.100.h5') as statepoint:
    photon_tally = statepoint.get_tally(name="photon_dose_on_mesh")

    # normalising this tally is a little different to other examples as the source strength has been using units of photons per second.
    # tally.mean is in units of pSv-cm3/source photon.
    # as source strength is in photons_per_second this changes units to pSv-/second

    # multiplication by pico_to_micro converts from (pico) pSv/s to (micro) uSv/s
    # dividing by mesh voxel volume cancles out the cm3 units
    # could do the mesh volume scaling on the plot and vtk functions but doing it here instead
    scaling_factor = (seconds_to_hours * pico_to_micro) / mesh.volumes[0][0][0]

    # plot = plot_mesh_tally(
    #         tally=photon_tally,
    #         basis="yz",
    #         # score='flux', # only one tally so can make use of default here
    #         value="mean",
    #         colorbar_kwargs={
    #             'label': "Decay photon dose [µSv/h]",
    #         },
    #         outl
    #         norm=LogNorm(),
    #         volume_normalization=False,  # this is done in the scaling_factor
    #         scaling_factor=scaling_factor,
    #     )
    plot = plot_mesh_tally(
        basis="xy",  # as the mesh dimention is [1,40,40] only the yz basis can be plotted
        tally=photon_tally,
        outline=True,  # enables an outline around the geometry
        geometry=model.geometry,  # needed for outline
        norm=LogNorm(),  # log scale
        colorbar=False,
    )
    plot.figure.savefig(f"activation_dose_map_timestep_{100}")




