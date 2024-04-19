import os
import openmc
import numpy as np

from .materials import dt_plasma, flibe, burner_mixture, v4cr4ti, tungsten
from barc_blanket.materials.waste_classification import separate_nuclides

# Default model parameters
# TODO: this all assumes a circular cross-section, which is not necessarily the case
# Must determine if this is a reasonable assumption
# See 'simple_toroidal.png' for a diagram of the geometry

DEFAULT_PARAMETERS = {
    'major_radius': 450,            # All dimensions are in cm
    'plasma_minor_radius': 138,
    'sol_width': 2,
    'first_wall_thickness': 0.1,          # How thick the plasma facing material is
    'vacuum_vessel_thickness': 1,         # How thick the vacuum vessel is
    'cooling_channel_width': 1,           # Width of the flowing coolant
    'cooling_vessel_thickness': 1,        # How thick the cooling vessel is
    'blanket_width': 130,                 # Width of the bulk molten salt blanket
    'blanket_vessel_thickness': 8,        # How thick the blanket vessel is

    'section_angle': 45,            # Angle of the toroidal section in degrees

    'li6_enrichment': 0.076,        # atom% enrichment of Li6 in the FLiBe
    'slurry_ratio': 0.01,           # weight% slurry in the burner blanket
    'removed_Sr90': 0.0,            # fraction of Sr90 removed from the waste stream
    'removed_Cs137': 0.0,           # fraction of Cs137 removed from the waste stream
    'removed_Tc99': 0.0,            # fraction of Tc99 removed from the waste stream
  
    'batches': 50,               # Number of batches to run
    'particles': int(1e5),       # Number of particles per batch

    'photon_transport': False    # Whether to run photon transport
}

def make_model(new_model_config=None):
    """Create an OpenMC model using the given configuration
    
    Parameters:
    ----------
    new_model_config : dict, optional
        Dictionary containing the model configuration.
        If not provided, the values listed in DEFAULT_PARAMETERS will be used.

    Returns:
    -------
    model : openmc.Model
        An OpenMC model object
    """

    if new_model_config is None:
        model_config = DEFAULT_PARAMETERS
    else:
        model_config = new_model_config.copy()
        for key in DEFAULT_PARAMETERS:
            if key not in new_model_config:
                model_config[key] = DEFAULT_PARAMETERS[key]
            else:
                print(f"Using set value for {key}:\t {model_config[key]}")

    ######################
    ## Assign Materials ##
    ######################

    plasma_material = dt_plasma()
    first_wall_material = tungsten()
    vacuum_vessel_material = v4cr4ti()
    flibe_material = flibe(model_config['li6_enrichment'])
    cooling_channel_material = burner_mixture(model_config['slurry_ratio'], flibe=flibe_material)
    cooling_vessel_material = v4cr4ti()
    blanket_material = burner_mixture(model_config['slurry_ratio'], flibe=flibe_material)
    blanket_vessel_material = v4cr4ti()

    # Remove some nuclides from blanket materials if applicable
    removed_materials_dict = {}
    if not np.isclose(model_config['removed_Sr90'], 0):
        removed_materials_dict['Sr90'] = model_config['removed_Sr90']
    if not np.isclose(model_config['removed_Cs137'], 0):
        removed_materials_dict['Cs137'] = model_config['removed_Cs137']
    if not np.isclose(model_config['removed_Tc99'], 0):
        removed_materials_dict['Tc99'] = model_config['removed_Tc99']

    if len(removed_materials_dict.keys()) > 0:
        cooling_channel_material = separate_nuclides(cooling_channel_material, removed_materials_dict)
        blanket_material = separate_nuclides(blanket_material, removed_materials_dict)

    #####################
    ## Define Geometry ##
    #####################

    R = model_config['major_radius']
    a = model_config['plasma_minor_radius']
    sol_width = model_config['sol_width']
    first_wall_thickness = model_config['first_wall_thickness']

    vacuum_vessel_thickness = model_config['vacuum_vessel_thickness']
    cooling_channel_width = model_config['cooling_channel_width']
    cooling_vessel_thickness = model_config['cooling_vessel_thickness']
    blanket_width = model_config['blanket_width']
    blanket_vessel_thickness = model_config['blanket_vessel_thickness']

    plasma_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=a,c=a)
    
    first_wall_inner_radius = a + sol_width
    first_wall_outer_radius = first_wall_inner_radius + first_wall_thickness
    first_wall_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=first_wall_inner_radius,c=first_wall_inner_radius)
    first_wall_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=first_wall_outer_radius,c=first_wall_outer_radius)

    vacuum_vessel_inner_radius = first_wall_outer_radius # These two are in contact
    vacuum_vessel_outer_radius = vacuum_vessel_inner_radius + vacuum_vessel_thickness
    vacuum_vessel_inner_surface = first_wall_outer_surface
    vacuum_vessel_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=vacuum_vessel_outer_radius,c=vacuum_vessel_outer_radius)

    cooling_vessel_inner_radius = vacuum_vessel_outer_radius+cooling_channel_width
    cooling_vessel_outer_radius = cooling_vessel_inner_radius + cooling_vessel_thickness
    cooling_vessel_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=cooling_vessel_inner_radius,c=cooling_vessel_inner_radius)
    cooling_vessel_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=cooling_vessel_outer_radius,c=cooling_vessel_outer_radius)

    blanket_vessel_inner_radius = cooling_vessel_outer_radius+blanket_width
    blanket_vessel_outer_radius = blanket_vessel_inner_radius + blanket_vessel_thickness
    blanket_vessel_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=blanket_vessel_inner_radius,c=blanket_vessel_inner_radius)
    blanket_vessel_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=blanket_vessel_outer_radius,c=blanket_vessel_outer_radius)

    bounding_sphere_surface = openmc.Sphere(r=2*R, boundary_type="vacuum")

    # Make two planes to cut the torus into a section
    # Angle follows right hand rule around z axis (https://www.desmos.com/3d/214a6bb908)
    section_angle_rad = np.radians(model_config['section_angle'])
    x_coeff, y_coeff = np.sin(section_angle_rad), -np.cos(section_angle_rad)
    xz_plane = openmc.Plane(a=0, b=1, boundary_type='periodic')
    angled_plane = openmc.Plane(a=x_coeff, b=y_coeff, boundary_type='periodic')
    xz_plane.periodic_surface = angled_plane
    torus_section = +xz_plane & +angled_plane

    volume_correction = model_config['section_angle']/360

    plasma_cell = openmc.Cell(
        name='plasma_cell',
        region=-plasma_surface & torus_section,
        fill=plasma_material
    )

    sol_cell = openmc.Cell(
        name='sol_cell',
        region=+plasma_surface & -first_wall_inner_surface & torus_section,
        fill=None
    )

    first_wall_cell = openmc.Cell(
        name='first_wall_cell',
        region=+first_wall_inner_surface & -vacuum_vessel_inner_surface & torus_section,
        fill=first_wall_material
    )
    first_wall_cell.fill.volume = (2*np.pi*R)*np.pi*(first_wall_outer_radius**2 - first_wall_inner_radius**2)*volume_correction

    vacuum_vessel_cell = openmc.Cell(
        name='vacuum_vessel_cell',
        region=+vacuum_vessel_inner_surface & -vacuum_vessel_outer_surface & torus_section,
        fill=vacuum_vessel_material
    )
    vacuum_vessel_cell.fill.volume = (2*np.pi*R)*np.pi*(vacuum_vessel_outer_radius**2 - vacuum_vessel_inner_radius**2)*volume_correction

    cooling_channel_cell = openmc.Cell(
        name='cooling_channel_cell',
        region=+vacuum_vessel_outer_surface & -cooling_vessel_inner_surface & torus_section,
        fill=cooling_channel_material
    )
    cooling_channel_cell.fill.volume = (2*np.pi*R)*np.pi*(cooling_vessel_inner_radius**2 - vacuum_vessel_outer_radius**2)*volume_correction

    cooling_vessel_cell = openmc.Cell(
        name='cooling_vessel_cell',
        region=+cooling_vessel_inner_surface & -cooling_vessel_outer_surface & torus_section,
        fill=cooling_vessel_material
    )
    cooling_vessel_cell.fill.volume = (2*np.pi*R)*np.pi*(cooling_vessel_outer_radius**2 - cooling_vessel_inner_radius**2)*volume_correction
    

    blanket_cell = openmc.Cell(
        name='blanket_cell',
        region=+cooling_vessel_outer_surface & -blanket_vessel_inner_surface & torus_section,
        fill=blanket_material
    )
    blanket_cell.fill.volume = (2*np.pi*R)*np.pi*(blanket_vessel_inner_radius**2 - cooling_vessel_outer_radius**2)*volume_correction

    blanket_vessel_cell = openmc.Cell(
        name='blanket_vessel_cell',
        region=+blanket_vessel_inner_surface & -blanket_vessel_outer_surface & torus_section,
        fill=blanket_vessel_material
    )
    blanket_vessel_cell.fill.volume = (2*np.pi*R)*np.pi*(blanket_vessel_outer_radius**2 - blanket_vessel_inner_radius**2)*volume_correction

    bounding_sphere_cell = openmc.Cell(
        name='bounding_sphere_cell',
        region=+blanket_vessel_outer_surface & -bounding_sphere_surface & torus_section,
        fill=None
    )

    universe = openmc.Universe()
    universe.add_cell(plasma_cell)
    universe.add_cell(sol_cell)
    universe.add_cell(first_wall_cell)
    universe.add_cell(vacuum_vessel_cell)
    universe.add_cell(cooling_channel_cell)
    universe.add_cell(cooling_vessel_cell)
    universe.add_cell(blanket_cell)
    universe.add_cell(blanket_vessel_cell)
    universe.add_cell(bounding_sphere_cell)
    geometry = openmc.Geometry(universe)

    #####################
    ## Define Settings ##
    #####################

    source = openmc.IndependentSource()
    source.particle = 'neutron'
    # This source shape is a thin wire in the plasma core
    radius = openmc.stats.Discrete([R], [1])
    z_values = openmc.stats.Discrete([0], [1])
    angle = openmc.stats.Uniform(a=np.radians(0), b=section_angle_rad)
    source.space = openmc.stats.CylindricalIndependent(
        r=radius, phi=angle, z=z_values, origin=(0., 0., 0.))
    source.angle = openmc.stats.Isotropic() # Isotropic directio neutron is launched
    source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

    settings = openmc.Settings(run_mode='fixed source')
    settings.photon_transport = model_config['photon_transport']
    settings.source = source
    settings.batches = model_config['batches']
    settings.particles = int(model_config['particles'])
    # Make statepoints every 5 batches, ensuring the final batch is always included
    statepoint_set = set([i for i in range(5, model_config['batches']+1, 5)])
    statepoint_set.add(model_config['batches'])
    settings.statepoint = {'batches': list(statepoint_set)}
    settings.output = {'tallies': True}

    #####################
    ## Define Tallies  ##
    #####################

    first_wall_cell_filter = openmc.CellFilter([first_wall_cell])
    vacuum_vessel_cell_filter = openmc.CellFilter([vacuum_vessel_cell])
    cooling_channel_cell_filter = openmc.CellFilter([cooling_channel_cell])
    cooling_vessel_cell_filter = openmc.CellFilter([cooling_vessel_cell])
    blanket_cell_filter = openmc.CellFilter([blanket_cell])
    tritium_cell_filter = openmc.CellFilter([cooling_channel_cell, blanket_cell])
    energy_filter = openmc.EnergyFilter(np.logspace(0,7)) # 1eV to 100MeV

    # mesh tally - flux
    tally1 = openmc.Tally(tally_id=1, name="flux_blanket")
    tally1.filters = [blanket_cell_filter,energy_filter]
    tally1.scores = ["flux"]

    # tbr
    tally2 = openmc.Tally(tally_id=2, name="tbr")
    tally2.filters = [tritium_cell_filter]
    tally2.scores = ["(n,Xt)"]

    #power deposition - heating
    # TODO: talk about doing coupled photon transport in the future
    # Definitely want the option to turn coupled photon on/off because it is computationally expensive
    first_wall_heating_tally = openmc.Tally(tally_id=3, name="neutron_heating_first_wall")
    first_wall_heating_tally.filters = [first_wall_cell_filter]
    first_wall_heating_tally.scores = ["heating"]

    vacuum_vessel_heating_tally = openmc.Tally(tally_id=4, name="neutron_heating_vacuum_vessel")
    vacuum_vessel_heating_tally.filters = [vacuum_vessel_cell_filter]
    vacuum_vessel_heating_tally.scores = ["heating"]

    cooling_channel_heating_tally = openmc.Tally(tally_id=5, name="neutron_heating_cooling_channel")
    cooling_channel_heating_tally.filters = [cooling_channel_cell_filter]
    cooling_channel_heating_tally.scores = ["heating"]

    cooling_vessel_heating_tally = openmc.Tally(tally_id=6, name="neutron_heating_cooling_vessel")
    cooling_vessel_heating_tally.filters = [cooling_vessel_cell_filter]
    cooling_vessel_heating_tally.scores = ["heating"]


    tallies = openmc.Tallies([tally1,tally2,
                              first_wall_heating_tally,vacuum_vessel_heating_tally,
                              cooling_channel_heating_tally,cooling_vessel_heating_tally]) 

    # Create model
    model = openmc.Model(geometry=geometry,
                         settings=settings, 
                         tallies=tallies)

    return model