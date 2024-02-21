# Simple file to ensure plotting cross sections works
import openmc
import barc_blanket.plot_cross_sections as pltxs

openmc.config['cross_sections'] = '/home/zkeith/proj/openmc_test/endfb-viii.0-hdf5/cross_sections.xml'

# 1. make a simple 'enriched uranium' material to compare results with
uranium_fuel = openmc.Material(name='Uranium Fuel')
uranium_fuel.add_nuclide('U235', 0.05)
uranium_fuel.add_nuclide('U238', 0.95)
uranium_fuel.add_element('O', 2)
uranium_fuel.set_density('g/cm3', 10.0)

# 2. Make a slurry of FLiBe and tank contents
flibe = openmc.Material(name="FLiBe")
flibe.add_element("Li", 2.0, "ao")
flibe.add_element("Be", 1.0, "ao")
flibe.add_element("F", 4.0, "ao")
flibe.set_density("g/cm3", 1.94)

# Just a guess until we have a better idea of what's in the tanks
tank_contents = openmc.Material(name='Tank Contents')
tank_contents.add_nuclide('Cs137', 0.2)
tank_contents.add_nuclide('Sr90', 0.2)
tank_contents.add_nuclide('Pu240', 0.2)
tank_contents.add_nuclide('Pu239', 0.2)
tank_contents.add_nuclide('Pu238', 0.2)
tank_contents.set_density('g/cm3', 8.0)

waste_fraction = 0.01 # Starting with 1% of the burner blanket contents being waste
slurry = openmc.Material(name="Slurry")
slurry.mix_materials([flibe, tank_contents], [1-waste_fraction, waste_fraction])
slurry.set_density("g/cm3", 4) # Just a guess, I think this is causing divide by zero error if you don't include it?


# 3. Plot the cross sections for the materials

fig = pltxs.nu_fission([uranium_fuel], normalize=True)

fig.savefig('nu_fission.png')
