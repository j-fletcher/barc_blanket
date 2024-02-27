from barc_blanket.vessel_activation import make_model, run_independent_vessel_activation, extract_activities
import os
import matplotlib.pyplot as plt

# Create a place to put all the files we'll be working with for depletion (TODO: standardize this structure with the rest of the group)
working_directory = "depletion"
if not os.path.exists(working_directory):
    os.makedirs(working_directory)
os.chdir(working_directory)

model = make_model()

if not os.path.exists("depletion_results.h5"):
    run_independent_vessel_activation(model, days=365*10, num_timesteps=50, source_rate=2e20)

timesteps, activities, dists = extract_activities()

# Plot the activities over time

fig, ax = plt.subplots()
ax.plot(timesteps, activities)
ax.set_xlabel("Time (days)")
ax.set_ylabel("Activity (Bq)")
ax.set_title("Vessel Activation")

# Save figure
fig.savefig("vessel_activation.png")