import os
import openmc
import pickle as pkl

from barc_blanket.models.barc_model_final import make_model, SECTION_CORRECTION
from barc_blanket.utilities import working_directory
from barc_blanket.materials.blanket_depletion import gw_to_neutron_rate

JOULES_PER_EV = 1.6e-19
FUSION_POWER_GW = 2.2

NUM_BATCHES = 50

def print_neutron_heating(tallies, cell_name_base, tally_id, model):

    source_particles_per_second = gw_to_neutron_rate(FUSION_POWER_GW, SECTION_CORRECTION)

    cell_name = f"{cell_name_base}_cell"
    midpl_name = f"{cell_name_base}_midpl"

    
    cell_volume_m3 = next(iter(model._cells_by_name[cell_name])).fill.volume / 1e6
    midpl_volume_m3 = next(iter(model._cells_by_name[midpl_name])).fill.volume / 1e6
    total_volume_m3 = cell_volume_m3 + midpl_volume_m3

    cell_neutron_heating_ev = tallies[tally_id].mean[0][0][0]
    midpl_neutron_heating_ev = tallies[tally_id+20].mean[0][0][0]

    cell_neutron_heating_joules = cell_neutron_heating_ev * JOULES_PER_EV / SECTION_CORRECTION
    midpl_neutron_heating_joules = midpl_neutron_heating_ev * JOULES_PER_EV / SECTION_CORRECTION

    cell_neutron_heating_watts = cell_neutron_heating_joules * source_particles_per_second
    midpl_neutron_heating_watts = midpl_neutron_heating_joules * source_particles_per_second

    total_neutron_heating_watts = cell_neutron_heating_watts + midpl_neutron_heating_watts

    midpl_neutron_heating_watts_per_m3 = midpl_neutron_heating_watts / midpl_volume_m3
    total_neutron_heating_watts_per_m3 = total_neutron_heating_watts / total_volume_m3

    peaking_factor = midpl_neutron_heating_watts_per_m3 / total_neutron_heating_watts_per_m3

    print(f"{cell_name_base} total neutron heating: {total_neutron_heating_watts/1e6} [MW]")
    print(f"{cell_name_base} volumetric neutron heating: {(total_neutron_heating_watts/1e6) / total_volume_m3} [MW/m3]")
    print(f"{cell_name_base} neutron heating peaking factor: {peaking_factor}")

    result_dict = {'total_MW': total_neutron_heating_watts/1e6,
                   'volumetric_MW_m3': (total_neutron_heating_watts/1e6) / total_volume_m3,
                   'peaking_factor': peaking_factor}
    
    return result_dict

with working_directory("neutron_heating_photon_transport"):

    model = make_model({
        'batches': NUM_BATCHES,
        'particles': 1e6,
        'photon_transport': True,
        'midplane_split': True
    })

    rerun_model = True
    if rerun_model is True:
        model.run()

    final_statepoint = openmc.StatePoint(f"statepoint.{NUM_BATCHES}.h5")

    # Get tally results
    # TODO: do this programmatically instead of hardcoded here
    tallies = final_statepoint.tallies

    all_results = {}

    all_results["first_wall"] = print_neutron_heating(tallies, "first_wall", 3, model)
    all_results["cooling_channel"] = print_neutron_heating(tallies, "cooling_channel", 4, model)
    all_results["cooling_vessel"] = print_neutron_heating(tallies, "cooling_vessel", 5, model)
    all_results["vacuum_vessel"] = print_neutron_heating(tallies, "vacuum_vessel", 6, model)
    all_results["blanket_vessel"] = print_neutron_heating(tallies, "blanket_vessel", 8, model)

    with open('neutron_heating_results.pkl', 'wb') as f:
        pkl.dump(all_results, f)