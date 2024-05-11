import os
import openmc

from barc_blanket.models.barc_model_final import make_model
from barc_blanket.utilities import working_directory
from barc_blanket.vessel_activation import get_cell_volume_in_mesh

JOULES_PER_EV = 1.6e-19

with working_directory("neutron_heating"):

    model = make_model({
        'batches': 10,
        'particles': 100
    })

    rerun_model = False
    if rerun_model is True:
        model.run()

    #get_cell_volume_in_mesh(model, n_samples=1000)

    final_statepoint = openmc.StatePoint("statepoint.20.h5")

    # Get tally results
    # TODO: do this programmatically instead of hardcoded here
    tallies = final_statepoint.tallies

    source_particles = 

    neutron_heating_first_wall_id = 3
    neutron_heating_first_wall_ev = tallies[neutron_heating_first_wall_id].mean[0][0][0]
    neutron_heating_first_wall_joules = neutron_heating_first_wall_ev * JOULES_PER_EV
    neutron_heating_first_wall_watts_per_gw = neutron_heating_first_wall_joules * SOURCE_PARTICLES_PER_GW
    print(f"First wall total neutron heating: {neutron_heating_first_wall_watts_per_gw/1e6} MW/GW")

    neutron_heating_vacuum_vessel_id = 4
    neutron_heating_vacuum_vessel_ev = tallies[neutron_heating_vacuum_vessel_id].mean[0][0][0]
    neutron_heating_vacuum_vessel_joules = neutron_heating_vacuum_vessel_ev * JOULES_PER_EV
    neutron_heating_vacuum_vessel_watts_per_gw = neutron_heating_vacuum_vessel_joules * SOURCE_PARTICLES_PER_GW
    print(f"Vacuum vessel neutron heating: {neutron_heating_vacuum_vessel_watts_per_gw/1e6} MW/GW")

    neutron_heating_cooling_channel_id = 5
    neutron_heating_cooling_channel_ev = tallies[neutron_heating_cooling_channel_id].mean[0][0][0]
    neutron_heating_cooling_channel_joules = neutron_heating_cooling_channel_ev * JOULES_PER_EV
    neutron_heating_cooling_channel_watts_per_gw = neutron_heating_cooling_channel_joules * SOURCE_PARTICLES_PER_GW
    print(f"Cooling channel neutron heating: {neutron_heating_cooling_channel_watts_per_gw/1e6} MW/GW")

    neutron_heating_cooling_vessel_id = 6
    neutron_heating_cooling_vessel_ev = tallies[neutron_heating_cooling_vessel_id].mean[0][0][0]
    neutron_heating_cooling_vessel_joules = neutron_heating_cooling_vessel_ev * JOULES_PER_EV
    neutron_heating_cooling_vessel_watts_per_gw = neutron_heating_cooling_vessel_joules * SOURCE_PARTICLES_PER_GW
    print(f"Cooling vessel neutron heating: {neutron_heating_cooling_vessel_watts_per_gw/1e6} MW/GW")

    
    # TODO do some volume weighting or whatever to get an actual TBR
    