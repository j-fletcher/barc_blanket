# %%
import openmc
import openmc_data_downloader as odd
import numpy as np

# this script uses openmc data downloader for reproducability 

# %%  Define Materials

# Plasma
dt_plasma = openmc.Material(name='dt_plasma')
dt_plasma.add_nuclide('H2', 1.0)
dt_plasma.add_nuclide('H3', 1.0)
dt_plasma.set_density('g/cm3', 1e-5)

# FLIBE
flibe = openmc.Material(name="flibe")
flibe.add_element("Li", 2.0, "ao")
flibe.add_element("Be", 1.0, "ao")
flibe.add_element("F", 4.0, "ao")
flibe.set_density("g/cm3", 1.94)

# Inconel 718 -
inconel718 = openmc.Material(name='inconel718')
inconel718.add_element('Ni', 53.0, 'wo')
inconel718.add_element('Cr', 19.06, 'wo')
inconel718.add_element('Nb', 5.08, 'wo')
inconel718.add_element('Mo', 3.04, 'wo')
inconel718.add_element('Ti', 0.93, 'wo')
inconel718.add_element('Al', 0.52, 'wo')
inconel718.add_element('Co', 0.11, 'wo')
inconel718.add_element('Cu', 0.02, 'wo')
inconel718.add_element('C', 0.021, 'wo')
inconel718.add_element('Fe', 18.15, 'wo')
inconel718.set_density('g/cm3', 8.19)

# Eurofer
eurofer = openmc.Material(name='eurofer')
eurofer.add_element('Cr', 8.99866, 'wo')
eurofer.add_element('C', 0.109997, 'wo')
eurofer.add_element('W', 1.5, 'wo')
eurofer.add_element('V', 0.2, 'wo')
eurofer.add_element('Ta', 0.07, 'wo')
eurofer.add_element('B', 0.001, 'wo')
eurofer.add_element('N', 0.03, 'wo')
eurofer.add_element('O', 0.01, 'wo')
eurofer.add_element('S', 0.001, 'wo')
eurofer.add_element('Fe', 88.661, 'wo')
eurofer.add_element('Mn', 0.4, 'wo')
eurofer.add_element('P', 0.005, 'wo')
eurofer.add_element('Ti', 0.01, 'wo')
eurofer.set_density('g/cm3', 7.798)

# V-4Cr-4Ti - pure -(from Segantin TRE https://github.com/SteSeg/tokamak_radiation_environment)
v4cr4ti = openmc.Material(name='v4cr4ti')
v4cr4ti.add_element('V', 0.92, 'wo')
v4cr4ti.add_element('Cr', 0.04, 'wo')
v4cr4ti.add_element('Ti', 0.04, 'wo')
v4cr4ti.set_density('g/cm3', 6.06)

#UF4 Uranium Tetrafuoride
uf4 = openmc.Material(name='uf4')
uf4.add_nuclide('U235',1.0)
uf4.add_element('F',4.0)
uf4.set_density('g/cm3',6.7)

#FLiBe UF4 Mixture
fueled_flibe = openmc.Material.mix_materials([flibe,uf4],[0.99,0.01],'vo')

# %% Geometry Definitions

# All measurements in cm
R = 680 # major radius
a = 120 # plasma minor radius
sol_thickness = 125 # distance between plasma boundary and 
vv_thickness = 127 # vaccuum vessel thickness
fusion_blanket_out = 142
burner_blanket_in = 147 # currently I have the two vessels separated by 5cm vanadium alloy
burner_blanket_out = 247

bounding_sphere = openmc.Sphere(r=2*R, boundary_type="vacuum")

plasma = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=a,c=a)

vv = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=sol_thickness,c=sol_thickness)

blanket_tank = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=vv_thickness,c=vv_thickness)

tank_divider_in = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=fusion_blanket_out,c=fusion_blanket_out)

tank_divider_out = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=burner_blanket_in,c=burner_blanket_in)

burner_tank = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=burner_blanket_out,c=burner_blanket_out)
# %% Cell definitions
plasma_cell = openmc.Cell()
plasma_cell.region = -plasma
plasma_cell.fill = dt_plasma

# this is just an empty cell separating plasma from vv
sol_cell = openmc.Cell()
sol_cell.region = +plasma & -vv

vv_cell = openmc.Cell()
vv_cell.region = +vv & -blanket_tank
vv_cell.fill = v4cr4ti #vanadium (switched from inconel718)

blanket1_cell = openmc.Cell()
blanket1_cell.region = +blanket_tank & -tank_divider_in
blanket1_cell.fill = flibe

tank_divider_cell = openmc.Cell()
tank_divider_cell.region = +tank_divider_in & -tank_divider_out
tank_divider_cell.fill = v4cr4ti

blanket2_cell = openmc.Cell()
blanket2_cell.region = +tank_divider_out & -burner_tank
blanket2_cell.fill = fueled_flibe

container_cell = openmc.Cell()
container_cell.region = +burner_tank & -bounding_sphere

universe = openmc.Universe()
universe.add_cell(plasma_cell)
universe.add_cell(sol_cell)
universe.add_cell(vv_cell)
universe.add_cell(blanket1_cell)
universe.add_cell(tank_divider_cell)
universe.add_cell(blanket2_cell)
universe.add_cell(container_cell)

# %% Source Definition

source = openmc.IndependentSource()
source.particle = 'neutron'
radius = openmc.stats.Discrete([680], [1]) # centered at major radius
z_values = openmc.stats.Discrete([0], [1])
angle = openmc.stats.Uniform(a=np.radians(0), b=np.radians(360))
source.space = openmc.stats.CylindricalIndependent(
    r=radius, phi=angle, z=z_values, origin=(0., 0., 0.))
source.angle = openmc.stats.Isotropic()
source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

# %% Settings

settings = openmc.Settings(run_mode='fixed source')
settings.photon_transport = False
settings.source = source
settings.batches = 50
settings.particles = int(1e5) # modify this to shorten simulation, default was 1e6 
settings.statepoint = {'batches': [
    5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]}
settings.output = {'tallies': True}

# %% Filters and Tallies

burner_cell_filter = openmc.CellFilter(blanket2_cell) # just burner blanket
tbr_cell_filter = openmc.CellFilter([blanket1_cell,blanket2_cell]) # includes both fusion and burner
energy_filter = openmc.EnergyFilter(np.logspace(0,7)) # 1eV to 100MeV

# mesh tally - flux

tally1 = openmc.Tally(tally_id=1, name="flux_burner")
tally1.filters = [burner_cell_filter,energy_filter]
tally1.scores = ["flux"]

# tbr

tally2 = openmc.Tally(tally_id=2, name="tbr")
tally2.filters = [tbr_cell_filter]
tally2.scores = ["(n,Xt)"]

#power deposition - heating-local

tally3 = openmc.Tally(tally_id=3, name="heating_burner")
tally3.filters = [burner_cell_filter]
tally3.scores = ["heating-local"]
# %% 

# Define materials, geomery, and tallies
materials = openmc.Materials([dt_plasma,flibe,inconel718,eurofer,v4cr4ti,fueled_flibe])
geometry = openmc.Geometry(universe)
tallies = openmc.Tallies([tally1,tally2,tally3])    

# Download cross cection data for materials
materials.download_cross_section_data(
        libraries=['ENDFB-7.1-NNDC', 'TENDL-2019'],
        set_OPENMC_CROSS_SECTIONS=True,
        particles=["neutron"],
    )

model = openmc.Model(geometry=geometry,
                     settings=settings, tallies=tallies)

#model = openmc.Model(materials=materials, geometry=geometry,
#                     settings=settings, tallies=tallies)

# Export and run model
model.export_to_model_xml()
model.run(threads=8, geometry_debug=True)