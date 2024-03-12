from barc_blanket.vessel_activation import run_independent_vessel_activation, extract_activities, extract_nuclides
import os
import matplotlib.pyplot as plt
from barc_blanket.utilities import working_directory
from barc_blanket.models.barc_model_simple_toroidal import make_model

# Create a place to put all the files we'll be working with for depletion
with working_directory("independent_vessel_decay"):

    model_config = {"batches": 10,
                    "particles": 1000,
                    "slurry_ratio": 0}
    model = make_model(model_config)
    # Save model as xml
    model.export_to_model_xml()

    rerun_depletion = True
    if not os.path.exists("depletion_results.h5") or rerun_depletion:
        run_independent_vessel_activation(model, days=365, num_timesteps=100)

    timesteps, nuclides = extract_nuclides(model)

    normalized_nuclides = {}
    for key in nuclides:
        normalized_nuclides[key] = (nuclides[key] - nuclides[key][0]) / nuclides[key][0]

    # Plot the change in nuclide concentration over time
    fig, ax = plt.subplots()
    for key in normalized_nuclides:
        ax.plot(timesteps, normalized_nuclides[key], label=key)

    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Normalized Concentration")
    ax.set_title("Nuclide Atom Count")
    ax.legend(loc='upper right')
    fig.savefig("nuclide_concentration.png")


    # # TODO: see how it's doing this
    timesteps, vv_activities = extract_activities(model, "vv_cell")
    timesteps, bv_activities = extract_activities(model, "bv_cell")

    # # Plot the activities over time

    fig, ax = plt.subplots()
    #ax.plot(timesteps, activities)
    ax.plot(timesteps, vv_activities[1], label="Vacuum Vessel")
    ax.plot(timesteps, bv_activities[1], label="Blanket Vessel")
    ax.legend()
    ax.set_xlabel("Time (days)")
    ax.set_ylabel("Activity (Bq)")
    ax.set_title("Vessel Activation")

    # Save figure
    fig.savefig("vessel_activation.png")

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