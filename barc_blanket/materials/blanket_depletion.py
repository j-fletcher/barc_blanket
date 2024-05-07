import numpy as np
import pickle as pkl
import matplotlib.pyplot as plt

import openmc.deplete
from barc_blanket.materials.waste_classification import sum_of_fractions, remove_flibe, remove_tritium
from barc_blanket.models.barc_model_final import SECTION_CORRECTION

def gw_to_neutron_rate(gw, section_correction=SECTION_CORRECTION):
    """Convert GW of fusion power to neutron rate in n/s
    
    Parameters
    ----------
    gw : float
        Fusion power in GW
    section_correction : float, optional
        Fraction of total torus that the section takes up to adjust total neutron rate.
        Default is using the section correction from barc_model_final.py
    """

    efus = 17.6e6  # eV
    ev2j = 1.60218e-19
    neutron_rate = gw*1e9 * section_correction / (efus * ev2j)  # n/s

    return neutron_rate

def run_coupled_depletion(model, timesteps_years, fusion_power):
    """ Run coupled depletion for a given model and timesteps
    Results are saved in 'depletion_results.h5' file in whatever directory called this function

    Parameters
    ----------
    model : openmc.model.Model
        Model to run depletion for
    timesteps_years : numpy.ndarray
        Array of timesteps to run depletion for (in years)
    fusion_power : float
        Fusion power in GW
    """

    timesteps_days = np.array(timesteps_years) * 365  # convert to days

    source_rates = np.ones_like(timesteps_years) * gw_to_neutron_rate(fusion_power)

    op = openmc.deplete.CoupledOperator(model, 
                                    reduce_chain=True, 
                                    reduce_chain_level=5, 
                                    normalization_mode='source-rate')
    
    openmc.deplete.CECMIntegrator(op, timesteps_days, source_rates=source_rates, timestep_units='d').integrate()

def postprocess_coupled_depletion(flibe_material_index, remove_C14=False):
    """Postprocess the results of a coupled depletion run
    
    Assumed to be ran in the same directory as the depletion results
    """

    # Load the results
    results = openmc.deplete.Results("depletion_results.h5")

    times_days = results.get_times()
    times_years = times_days / 365
    # round to nearest int
    times_years = np.round(times_years).astype(int)

    blanket_composition_at_time = []

    for i, time in enumerate(times_years):
        materials = results.export_to_materials(burnup_index=i, path='materials.xml')
        blanket_composition_at_time.append(materials[flibe_material_index])

    blanket_result_dictionary = {}
    for blanket_material, time in zip(blanket_composition_at_time, times_years):

        removed_tritium = remove_tritium(blanket_material, 0.9)
        removed_flibe = remove_flibe(removed_tritium, 0.9)
        sample_material = removed_flibe


        table_1_sum_of_fractions, table_1_culprits = sum_of_fractions(sample_material, 1, None, remove_C14=remove_C14)
        table_2_sum_of_fractions, table_2_culprits = sum_of_fractions(sample_material, 2, 3)

        print(f"Time: {time} years")
        print(f"Table 1 sum of fractions: {table_1_sum_of_fractions:0.2f}")
        print(f"Table 2 sum of fractions: {table_2_sum_of_fractions:0.2f}")

        blanket_result_dictionary[time] = {'table_1_sum_of_fractions': table_1_sum_of_fractions,
                                    'table_1_culprits': table_1_culprits,
                                    'table_2_sum_of_fractions': table_2_sum_of_fractions,
                                    'table_2_culprits': table_2_culprits}
        
    full_result_dictionary = {'blanket': blanket_result_dictionary}

    # Pickle the results
    with open('waste_classification_results.pkl', 'wb') as f:
        pkl.dump(full_result_dictionary, f)

def plot_results(case:str, print_name:str):

    with open(f'waste_classification_results.pkl', 'rb') as f:
        result_dictionary = pkl.load(f)

    for cell in result_dictionary.keys():
        print(f"Cell: {cell}")
        cell_result_dictionary = result_dictionary[cell]

        # Sum of fractions plot

        table_1_sum_of_fractions = []
        table_2_sum_of_fractions = []

        times = list(cell_result_dictionary.keys())

        for time in times:
            table_1_sum_of_fractions.append(cell_result_dictionary[time]['table_1_sum_of_fractions'])
            table_2_sum_of_fractions.append(cell_result_dictionary[time]['table_2_sum_of_fractions'])

        # Plot the results on a single plot with semilog y axis
        fig, ax = plt.subplots()
        ax.semilogy(times, table_1_sum_of_fractions, label='Table 1')
        ax.semilogy(times, table_2_sum_of_fractions, label='Table 2')
        ax.set_xlabel('Time (years)', fontsize=16)
        ax.set_ylabel('Sum of Fractions', fontsize=16)

        # Put a horizontal line at 1 for reference
        ax.axhline(1, color='black', linestyle='--', label='CCLLW', linewidth=2)
        # Put a horizontal line at 2.33 for reference
        ax.axhline(2.33, color='purple', linestyle='--', label='CCLLW with Vitrification', linewidth=2)
        ax.set_xlim(0, 100)
        ax.legend()
        ax.set_title(f'{print_name} Sum of Fractions', fontsize=18)

        # Save figure to file
        fig.savefig(f'{case}_sum_of_fractions.png')

        # Culprits in Table 1

        # For the last time step, plot the 6 largest contributors to the sum of fractions
        time = times[-1]
        table_1_culprits = cell_result_dictionary[time]['table_1_culprits']

        # Find the 6 largest entries
        table_1_culprits = sorted(table_1_culprits.items(), key=lambda x: x[1], reverse=True)[:6]

        fig, ax = plt.subplots()
        ax.bar([x[0] for x in table_1_culprits], [x[1] for x in table_1_culprits], label='Table 1')
        ax.set_xlabel('Nuclide', fontsize=16)
        ax.set_ylabel('Fraction', fontsize=16)

        ax.set_title(f'{print_name} Table 1 dominant nuclides at {time} years', fontsize=18)

        # Save figure to file
        fig.savefig(f'{case}_table_1_culprits.png')
        