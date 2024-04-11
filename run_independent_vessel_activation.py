from barc_blanket.vessel_activation import run_independent_vessel_activation, extract_activities, extract_nuclides
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from barc_blanket.utilities import working_directory
from barc_blanket.models.barc_model_simple_toroidal import make_model

# Create a place to put all the files we'll be working with for depletion
with working_directory("independent_vessel_activation"):

    model_config = {"batches": 10,
                    "particles": 1000,
                    "slurry_ratio": 0,
                    "section_angle": 10,}
    model = make_model(model_config)
    # Save model as xml
    model.export_to_model_xml()

    rerun_depletion = True
    times = np.geomspace(0.01, 365, 100)
    if not os.path.exists("depletion_results.h5") or rerun_depletion:
        run_independent_vessel_activation(model, times=times)

    nuclide_times, nuclides = extract_nuclides(model, cell_name="blanket_vessel_cell", nuclide_names=["V49"])

    # Plot the change in nuclide concentration over time
    fig, ax = plt.subplots()
    for key in nuclides:
        ax.plot(nuclide_times, nuclides[key], label=key)

    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Nuclide Count")
    ax.set_title("Nuclide Atom Count over Time")
    ax.legend(loc='upper right')
    fig.savefig("nuclide_concentration.png")

    activity_times, blanket_vessel_activities = extract_activities(model, "blanket_vessel_cell")

    # # Plot the activities over time

    plt.figure(figsize=(8, 6))
    plt.style.use('seaborn-v0_8-poster')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.grid(True, color='w', linestyle='-', linewidth=1.5)
    plt.gca().patch.set_facecolor('0.92')

    plt.plot(activity_times, blanket_vessel_activities[1], label="Blanket Vessel")
    plt.legend(loc='right')
    plt.xlabel("Time [days]")
    plt.ylabel("Activity [Bq]")
    plt.title("Vessel Activation")
    plt.yscale('symlog', linthresh=1e12)
    plt.xlim(0, max(activity_times))
    plt.ylim(0, max(blanket_vessel_activities[1])*1.01)
    

    # Save figure
    plt.savefig("blanket_vessel_activation.png")

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