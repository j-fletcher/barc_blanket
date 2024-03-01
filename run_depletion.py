from barc_blanket.vessel_activation import make_model, run_independent_vessel_activation, extract_activities, extract_nuclides
import os
import matplotlib.pyplot as plt

# Create a place to put all the files we'll be working with for depletion (TODO: standardize this structure with the rest of the group)
working_directory = "depletion"
if not os.path.exists(working_directory):
    os.makedirs(working_directory)
os.chdir(working_directory)

model = make_model()

rerun_depletion = True
if not os.path.exists("depletion_results.h5") or rerun_depletion:
    run_independent_vessel_activation(model, days=100, num_timesteps=10, source_rate=2e20)


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
# timesteps, activities = extract_activities(model)

# # Plot the activities over time

# fig, ax = plt.subplots()
# #ax.plot(timesteps, activities)
# ax.plot(timesteps, activities[1])
# ax.set_xlabel("Time (days)")
# ax.set_ylabel("Activity (Bq)")
# ax.set_title("Vessel Activation")

# # Save figure
# fig.savefig("vessel_activation.png")

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