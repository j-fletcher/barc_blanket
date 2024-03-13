from barc_blanket.vessel_activation import run_independent_vessel_activation, extract_activities, extract_nuclides
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from barc_blanket.utilities import working_directory
from barc_blanket.models.barc_model_simple_toroidal import make_model

# Create a place to put all the files we'll be working with for depletion
with working_directory("independent_vessel_activation"):

    model_config = {"batches": 10,
                    "particles": 1000,
                    "slurry_ratio": 0}
    model = make_model(model_config)
    # Save model as xml
    model.export_to_model_xml()

    rerun_depletion = False
    if not os.path.exists("depletion_results.h5") or rerun_depletion:
        run_independent_vessel_activation(model, days=365, num_timesteps=100)

    # timesteps, nuclides = extract_nuclides(model)

    # normalized_nuclides = {}
    # for key in nuclides:
    #     normalized_nuclides[key] = (nuclides[key] - nuclides[key][0]) / nuclides[key][0]

    # # Plot the change in nuclide concentration over time
    # fig, ax = plt.subplots()
    # for key in normalized_nuclides:
    #     ax.plot(timesteps, normalized_nuclides[key], label=key)

    # ax.set_xlabel("Time (days)")
    # ax.set_ylabel("Normalized Concentration")
    # ax.set_title("Nuclide Atom Count")
    # ax.legend(loc='upper right')
    # fig.savefig("nuclide_concentration.png")

    timesteps, first_wall_activities = extract_activities(model, "first_wall_cell")
    timesteps, vv_activities = extract_activities(model, "vv_cell")
    timesteps, bv_activities = extract_activities(model, "bv_cell")

    # # Plot the activities over time

    plt.figure(figsize=(8, 6))
    plt.style.use('seaborn-v0_8-poster')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(True, color='w', linestyle='-', linewidth=1.5)
    plt.gca().patch.set_facecolor('0.92')

    #ax.plot(timesteps, activities)
    plt.plot(timesteps, first_wall_activities[1], label="First Wall")
    plt.plot(timesteps, vv_activities[1], label="Vacuum Vessel")
    plt.plot(timesteps, bv_activities[1], label="Blanket Vessel")
    plt.legend(loc='right')
    plt.xlabel("Time [days]")
    plt.ylabel("Activity [Bq]")
    plt.title("Vessel Activation")
    plt.yscale('symlog', linthresh=1e12)
    plt.xlim(0, max(timesteps))
    plt.ylim(0, max(max(first_wall_activities[1]), max(vv_activities[1]), max(bv_activities[1]))*1.01)
    

    # Save figure
    plt.savefig("vessel_activation.png")

    plt.rcParams.update(mpl.rcParamsDefault)
    

    # # fig, ax = plt.subplots()
    # # for i in [0, 20, len(timesteps)-1]:
    # #     energies = dists[i].x
    # #     activities = dists[i].p
    # #     ax.scatter(energies, activities, label=f"Step {i}")
    # #     ax.set_xlabel("Energy (MeV)")
    # #     ax.set_ylabel("Activity (Bq)")

    # ax.legend()
    # ax.set_title("Energy Spectrum")
    # fig.savefig("energy_spectrum.png")