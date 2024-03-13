import openmc
import openmc.model
import openmc.deplete
from barc_blanket.vessel_activation import run_independent_vessel_decay, extract_activities, extract_nuclides
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from barc_blanket.utilities import working_directory
from barc_blanket.models.materials import water

result_directory = "independent_vessel_activation"

# Create a place to put all the files we'll be working with for depletion
with working_directory("independent_vessel_decay"):

    # Load model
    activated_model = openmc.model.Model.from_model_xml(f"../{result_directory}/model.xml")

    # Replace the blanket material with water
    blanket_cell = next(iter(activated_model._cells_by_name["blanket_cell"]))
    blanket_cell.fill = water()

    results = openmc.deplete.Results(f"../{result_directory}/depletion_results.h5")

    # Replace the first wall material with activated first wall material
    first_wall_cell = next(iter(activated_model._cells_by_name["first_wall_cell"]))
    # Get the activated first wall material
    #activated_first_wall_material = results[99].get_material(str(3))
    #first_wall_cell.fill = activated_first_wall_material

    # Replace the vacuum vessel material with activated vacuum vessel material
    vv_cell = next(iter(activated_model._cells_by_name["vv_cell"]))
    # Get the activated vacuum vessel material
    #activated_vv_material = results[99].get_material(str(4))
    #vv_cell.fill = activated_vv_material

    # Replace the blanket vessel material with activated blanket vessel material
    bv_cell = next(iter(activated_model._cells_by_name["bv_cell"]))
    # Get the activated blanket vessel material
    #activated_bv_material = results[99].get_material(str(8))
    #bv_cell.fill = activated_bv_material

    decay_model = openmc.model.Model(geometry=activated_model.geometry, 
                                     settings=activated_model.settings)
    
    rerun_depletion = False
    if not os.path.exists("depletion_results.h5") or rerun_depletion:
        #run_independent_vessel_activation(decay_model, days=365, num_timesteps=100, source_rate=0)
        run_independent_vessel_decay(decay_model, results, days=365, num_timesteps=100)

    # TODO: see how it's doing this
    timesteps, first_wall_activities = extract_activities(decay_model, "first_wall_cell")
    timesteps, vv_activities = extract_activities(decay_model, "vv_cell")
    timesteps, bv_activities = extract_activities(decay_model, "bv_cell")

    # # Plot the activities over time

    plt.figure(figsize=(8, 6))
    plt.style.use('seaborn-v0_8-poster')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(True, color='w', linestyle='-', linewidth=1.5)
    plt.gca().patch.set_facecolor('0.92')

    #ax.plot(timesteps, activities)
    plt.plot(timesteps-365, first_wall_activities[1], label="First Wall")
    plt.plot(timesteps-365, vv_activities[1], label="Vacuum Vessel")
    plt.plot(timesteps-365, bv_activities[1], label="Blanket Vessel")
    plt.legend(loc='upper right')
    plt.xlabel("Time [days]")
    plt.ylabel("Activity [Bq]")
    plt.title("Vessel Decay")
    plt.yscale('symlog', linthresh=1e12)
    plt.xlim(0, 365)
    #plt.ylim(0, max(max(first_wall_activities[1]), max(vv_activities[1]), max(bv_activities[1]))*1.01)
    

    # Save figure
    plt.savefig("vessel_activation.png")

    plt.rcParams.update(mpl.rcParamsDefault)

