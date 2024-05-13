import openmc
import openmc.model
import openmc.deplete
from barc_blanket.vessel_activation import run_independent_vessel_decay, extract_activities, extract_decay_heat
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from barc_blanket.utilities import working_directory
from barc_blanket.models.materials import water
from barc_blanket.models.barc_model_final import SECTION_CORRECTION

result_directory = "independent_vessel_activation"

def plot_activities(decay_model):
    # Plot the activities over time

    timesteps, first_wall_activities = extract_activities(decay_model, "first_wall_cell")
    timesteps, vacuum_vessel_activities = extract_activities(decay_model, "vacuum_vessel_cell")
    timesteps, blanket_vessel_activities = extract_activities(decay_model, "blanket_vessel_cell")

    plt.figure(figsize=(8, 6))
    plt.style.use('seaborn-v0_8-poster')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(True, color='w', linestyle='-', linewidth=1.5)
    plt.gca().patch.set_facecolor('0.92')

    #ax.plot(timesteps, activities)
    plt.plot(timesteps-365, first_wall_activities[1], label="First Wall")
    plt.plot(timesteps-365, vacuum_vessel_activities[1], label="Vacuum Vessel")
    plt.plot(timesteps-365, blanket_vessel_activities[1], label="Blanket Vessel")
    plt.legend(loc='upper right')
    plt.xlabel("Time [days]")
    plt.ylabel("Activity [Bq]")
    plt.title("Vessel Decay")
    plt.yscale('symlog', linthresh=1e12)
    plt.xlim(0, 365)
    #plt.ylim(0, max(max(first_wall_activities[1]), max(vacuum_vessel_activities[1]), max(blanket_vessel_activities[1]))*1.01)
    
    # Save figure
    plt.savefig("vessel_activation.png")

    plt.rcParams.update(mpl.rcParamsDefault)


# Create a place to put all the files we'll be working with for depletion
with working_directory("independent_vessel_decay"):

    # Load model
    activated_model = openmc.model.Model.from_model_xml(f"../{result_directory}/model.xml")

    # Replace the flibe materials with water
    blanket_cell = next(iter(activated_model._cells_by_name["blanket_cell"]))
    blanket_cell.fill = water()
    cooling_channel_cell = next(iter(activated_model._cells_by_name["cooling_channel_cell"]))
    cooling_channel_cell.fill = water()

    results = openmc.deplete.Results(f"../{result_directory}/depletion_results.h5")

    # Replace the first wall material with activated first wall material
    first_wall_cell = next(iter(activated_model._cells_by_name["first_wall_cell"]))
    # Get the activated first wall material
    #activated_first_wall_material = results[99].get_material(str(3))
    #first_wall_cell.fill = activated_first_wall_material

    cooling_channel_cell = next(iter(activated_model._cells_by_name["cooling_channel_cell"]))

    # Replace the vacuum vessel material with activated vacuum vessel material
    vacuum_vessel_cell = next(iter(activated_model._cells_by_name["vacuum_vessel_cell"]))
    # Get the activated vacuum vessel material
    #activated_vacuum_vessel_material = results[99].get_material(str(4))
    #vacuum_vessel_cell.fill = activated_vacuum_vessel_material

    # Replace the blanket vessel material with activated blanket vessel material
    blanket_vessel_cell = next(iter(activated_model._cells_by_name["blanket_vessel_cell"]))
    # Get the activated blanket vessel material
    #activated_blanket_vessel_material = results[99].get_material(str(8))
    #blanket_vessel_cell.fill = activated_blanket_vessel_material

    decay_model = openmc.model.Model(geometry=activated_model.geometry, 
                                     settings=activated_model.settings)
    
    rerun_depletion = False
    times = np.geomspace(0.01, 5, 100)
    #times = [0, 0.25, 0.5, 0.75, 1, 1.5, 2, 2.5, 3]
    if not os.path.exists("depletion_results.h5") or rerun_depletion:
        run_independent_vessel_decay(decay_model, results, times=times)

    heat_times, first_wall_decay_heat = extract_decay_heat(decay_model, "first_wall_cell")
    heat_times, cooling_vessel_decay_heat = extract_decay_heat(decay_model, "cooling_vessel_cell")
    heat_times, vacuum_vessel_decay_heat = extract_decay_heat(decay_model, "vacuum_vessel_cell")
    heat_times, blanket_vessel_decay_heat = extract_decay_heat(decay_model, "blanket_vessel_cell")

    # Convert times to days
    heat_times_days = np.array(heat_times) / (60*60*24)

    first_wall_decay_heat_total_MW = (np.array(first_wall_decay_heat)/1e6) / SECTION_CORRECTION
    cooling_vessel_decay_heat_total_MW = (np.array(cooling_vessel_decay_heat)/1e6) / SECTION_CORRECTION
    vacuum_vessel_decay_heat_total_MW = (np.array(vacuum_vessel_decay_heat)/1e6) / SECTION_CORRECTION
    blanket_vessel_decay_heat_total_MW = (np.array(blanket_vessel_decay_heat)/1e6) / SECTION_CORRECTION

    # Make a dataframe of the decay heat
    df = pd.DataFrame({"Time [days]": heat_times_days, 
                       "First Wall [MW]": first_wall_decay_heat_total_MW,
                       "Cooling Vessel [MW]": cooling_vessel_decay_heat_total_MW,
                       "Vacuum Vessel [MW]": vacuum_vessel_decay_heat_total_MW, 
                       "Blanket Vessel [MW]": blanket_vessel_decay_heat_total_MW, 
                       })
    df.to_csv("vessel_decay_heat.csv", index=False)

    plt.rcParams.update(mpl.rcParamsDefault)
    # # Plot the decay heat over time
    plt.figure(figsize=(8, 6))
    plt.style.use('seaborn-v0_8-poster')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(True, color='w', linestyle='-', linewidth=1.5)
    plt.gca().patch.set_facecolor('0.92')

    plt.plot(heat_times_days[1:], first_wall_decay_heat_total_MW[1:], label="First Wall", linewidth=4)
    plt.plot(heat_times_days[1:], cooling_vessel_decay_heat_total_MW[1:], label="Cooling Vessel", linewidth=4)
    plt.plot(heat_times_days[1:], vacuum_vessel_decay_heat_total_MW[1:], label="Vacuum Vessel", linewidth=4)
    #plt.plot(heat_times_days[1:], blanket_vessel_decay_heat_total_MW[1:], label="Blanket Vessel", width=2)
    plt.xlim(0, max(heat_times_days))
    plt.legend(loc='upper right')
    plt.xlabel("Time [days]", fontsize=24)
    plt.ylabel("Decay Heat [MW]", fontsize=24)
    plt.title("Decay Heating in Replaceable Components", fontsize=26)
    plt.ylim(0, 3)
    # Make x ticks bigger
    plt.xticks(fontsize=20)
    # Make y ticks bigger
    plt.yticks(fontsize=20)
    # Make plot layout tight
    plt.tight_layout()

    plt.savefig("vessel_decay_heat.png")

    plt.rcParams.update(mpl.rcParamsDefault)
    

