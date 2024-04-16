import os
import openmc

from barc_blanket.models.barc_model_simple_toroidal import make_model, DEFAULT_PARAMETERS
from barc_blanket.models.barc_model_tungsten_cooling_channel import make_model_tungsten_cooling
from barc_blanket.utilities import working_directory

JOULES_PER_EV = 1.60218e-19
SOURCE_PARTICLES_PER_GW = 3.55e20

with working_directory("neutron_heating"):

    # Make and run the model
    # model = make_model({"batches": 20,
    #                     "particles": 1e4, 
    #                     "photon_transport": True,
    #                     "slurry_ratio": 0, 
    #                     "section_angle": 10, 
    #                     "first_wall_thickness": 1,
    #                     })
    
    model = make_model_tungsten_cooling({"batches": 20,
                                         "particles": 1e4, 
                                        "photon_transport": True,
                                        "slurry_ratio": 0, 
                                        "section_angle": 10, 
                                        "first_wall_thickness": 0.3,
                                        "cooling_channel_width": 1,
                                        "cooling_vessel_thickness": 0.3,
                                        })

    model.run()
    final_statepoint = openmc.StatePoint("statepoint.20.h5")

    # Get tally results
    # TODO: do this programmatically instead of hardcoded here
    tallies = final_statepoint.tallies

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
    