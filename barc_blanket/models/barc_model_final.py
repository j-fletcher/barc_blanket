import copy
import openmc
import numpy as np
from .materials import dt_plasma, tungsten, v4cr4ti, flibe, ss316L, shield, magnetmat

DEFAULT_PARAMETERS = {

    'major_radius': 480,
    'minor_radius': 130,
    'elongation': 1.8,

    'first_wall_thickness': 0.3,
    'cooling_channel_width': 1.2,
    'cooling_vessel_thickness': 0.3,
    'vacuum_vessel_thickness': 1,
    'blanket_vessel_thickness': 8,

    'blanket_outboard_gap': 130,
    'blanket_inboard_gap': 100,

    'shield_thickness': 35,
    'thermal_shield_thickness': 3,

    'magnetcase_thickness': 5,
    'magnet_thickness': 65,

    'plasma_material': dt_plasma(),
    'first_wall_material': tungsten(),
    'cooling_vessel_material': tungsten(),
    'vacuum_vessel_material': v4cr4ti(),
    'blanket_vessel_material': v4cr4ti(),
    'blanket_material': flibe(),
    'shield_material': shield(),
    'magnetcase_material': ss316L(),
    'magnet_material': magnetmat(),

    'batches': 100,
    'inactive_batches': 0,
    'particles': int(1e6),

    'photon_transport': False,

    'midplane_split': False
}

BLANKET_MATERIAL_ID = 5
SECTION_CORRECTION = 1/14 # This model is a section of the torus, so the volume is 1/14 of the total volume

MIDPLANE_OFFSET = 20

def peaking_sector_volume(majorrad,c1,c2,b1,b2,z0=MIDPLANE_OFFSET,theta=2*np.pi*SECTION_CORRECTION):
    """Calculate the volume of the inboard midplane sector from a particluar layer, where the flux and heating peak.
       This is represented by the shaded region below.
  
       |          ## ##               
       |        #  # #  #             ^
       |       # #  |   # #        ^  |
       |     # #    |    # #       |  |
       |    # #     |     # #      b1 b2
    z0-|----%%%-----|-----# #------|  |
       |   %%%      |      # #     |  |
       |   %%%------0      # #     0  0
       |   %%%             # #
    -z0|----%%%-----------# #------
       |    # #           # #
       |     # #         # #
       |       # #      # #
       |        #  # #  #
       |          ## ## 
       |
       |      <--c1-0
       |   <-----c2-0
       |<--majorrad-0
        


    Parameters:
    ----------
    majorrad : float
        Toroidal major radius of the component [cm]
        i.e., R or blanket_vessel_major_radius
    c1 : float
        Semiminor axis of the inner ellipsoid poloidal surface [cm]
    c2 : float
        Semiminor axis of the outer ellipsoid poloidal surface [cm]
    b1 : float
        Semimajor axis of the inner ellipsoid poloidal surface [cm]
    b2 : float
        Semimajor axis of the outer ellipsoid poloidal surface [cm]
    z0 : float, optional
        Half-height [cm] of "midplane" sector. Larger z0 produces better tally statistics, but poorer fidelity
        when calculating peaking factors. Default = 10 cm
    theta : float, optional
        Angle [rad] of toroidal section. Default = 2pi/(number of sections in torus)

    Returns:
    -------

    sector_volume : float
        Volume of inboard midplane sector.
    """

    term1 = z0*(c1**2) - (z0**3)*(c1**2)/(3*b1**2)
    term2 = z0*(c2**2) - (z0**3)*(c2**2)/(3*b2**2)

    term3a = z0*np.sqrt((c1**2) - (c1**2)*(z0**2)/(b1**2))
    term3b1 = (b1*z0*np.sqrt((c1**2) - (c1**2)*(z0**2)/(b1**2)))/(c1*((b1**2)-(z0**2)))
    term3b = c1*b1*np.arctan(term3b1)
    term3 = majorrad*(term3a + term3b)

    term4a = z0*np.sqrt((c2**2) - (c2**2)*(z0**2)/(b2**2))
    term4b1 = (b2*z0*np.sqrt((c2**2) - (c2**2)*(z0**2)/(b2**2)))/(c2*((b2**2)-(z0**2)))
    term4b = c2*b2*np.arctan(term4b1)
    term4 = majorrad*(term4a + term4b)

    sector_volume = theta*(term1 - term2 - term3 + term4)

    return sector_volume

def total_layer_volume(majorrad,c1,c2,elongation,theta=2*np.pi*SECTION_CORRECTION):
    """Calculate the volume of a toroidal sector from a particluar layer.
       This returns the entire volume of the thick, elliptical toroidal shell section.
  
       |          ## ##               
       |        #  # #  #
       |       # #  |   # #
       |     # #    |    # #
       |    # #     |     # #
    z0-|----# #-----|-----# #------
       |   # #      |      # #
       |   # #------0      # #
       |   # #             # #
    -z0|----# #-----------# #------
       |    # #           # #
       |     # #         # #
       |       # #      # #
       |        #  # #  #
       |          ## ## 
       |
       |      <--c1-0
       |   <-----c2-0
       |<--majorrad-0
        


    Parameters:
    ----------
    majorrad : float
        Toroidal major radius of the component [cm]
        i.e., R or blanket_vessel_major_radius
    c1 : float
        Semiminor axis of the inner ellipsoid poloidal surface [cm]
    c2 : float
        Semiminor axis of the outer ellipsoid poloidal surface [cm]
    elongation : float
        Elliptical elongation factor
    theta : float, optional
        Angle [rad] of toroidal section. Default = 2pi/(number of sections in torus)
    Returns:
    -------

    layer_volume : float
        Volume of tokamak layer.
    """
    layer_volume = theta*majorrad*np.pi*(c2**2 - c1**2)*elongation

    return layer_volume

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

    #####################
    ## Define Geometry ##
    #####################

    R = model_config['major_radius']
    a = model_config['minor_radius']
    elongation = model_config['elongation']
    b=elongation*a

    first_wall_thickness = model_config['first_wall_thickness']
    cooling_channel_width = model_config['cooling_channel_width']
    cooling_vessel_thickness = model_config['cooling_vessel_thickness']
    vacuum_vessel_thickness = model_config['vacuum_vessel_thickness']
    blanket_vessel_thickness = model_config['blanket_vessel_thickness']
    neutron_shield_thickness = model_config['shield_thickness']
    inner_magnet_case_thickness = model_config['thermal_shield_thickness'] + model_config['magnetcase_thickness']
    # ^ this includes the thermal shield, being of the same material
    outer_magnet_case_thickness = model_config['magnetcase_thickness']
    magnet_coil_thickness = model_config['magnet_thickness']

    # Taking the provided minor radius as the smaller one
    # This is gonna make the thicknesses on the top and bottom a little bigger,
    # but I don't think that will make a significant difference
    plasma_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=b,c=a)

    first_wall_inner_radius = a
    first_wall_outer_radius = a + first_wall_thickness
    first_wall_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=first_wall_inner_radius*elongation,c=first_wall_inner_radius)
    first_wall_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=first_wall_outer_radius*elongation,c=first_wall_outer_radius)

    cooling_vessel_inner_radius = first_wall_outer_radius + cooling_channel_width
    cooling_vessel_outer_radius = cooling_vessel_inner_radius + cooling_vessel_thickness
    cooling_vessel_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=cooling_vessel_inner_radius*elongation,c=cooling_vessel_inner_radius)
    cooling_vessel_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=cooling_vessel_outer_radius*elongation,c=cooling_vessel_outer_radius)

    vacuum_vessel_inner_radius = cooling_vessel_outer_radius # These two are in contact
    vacuum_vessel_outer_radius = vacuum_vessel_inner_radius + vacuum_vessel_thickness
    vacuum_vessel_inner_surface = cooling_vessel_outer_surface
    vacuum_vessel_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=vacuum_vessel_outer_radius*elongation,c=vacuum_vessel_outer_radius)

    # Since the blanket is sorta offset by the gaps, I'm just modeling it as torus with a different major radius
    blanket_vessel_offset = (model_config['blanket_outboard_gap'] - model_config['blanket_inboard_gap'])/2
    blanket_vessel_average_width = (model_config['blanket_outboard_gap'] + model_config['blanket_inboard_gap'])/2
    blanket_vessel_major_radius = R + blanket_vessel_offset
    blanket_vessel_inner_radius = vacuum_vessel_outer_radius + blanket_vessel_average_width
    blanket_vessel_outer_radius = blanket_vessel_inner_radius + blanket_vessel_thickness
    blanket_vessel_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=blanket_vessel_major_radius,b=blanket_vessel_inner_radius*elongation*0.8,c=blanket_vessel_inner_radius)
    blanket_vessel_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=blanket_vessel_major_radius,b=blanket_vessel_outer_radius*elongation*0.8,c=blanket_vessel_outer_radius)

    # outer layers: neutron shield, thermal shield, magnet case, magnet coil, magnet case
    # the inner surface of the neutron shield is coincident with the blanket vessel
    neutron_shield_outer_radius = blanket_vessel_outer_radius + neutron_shield_thickness
    neutron_shield_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=blanket_vessel_major_radius,b=neutron_shield_outer_radius*elongation*0.8,c=neutron_shield_outer_radius)
    # in this case, the thermal shield and magnet case are also combined because they are of the same material
    inner_case_outer_radius = neutron_shield_outer_radius + inner_magnet_case_thickness
    inner_case_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=blanket_vessel_major_radius,b=inner_case_outer_radius*elongation*0.8,c=inner_case_outer_radius)
    magnet_coil_outer_radius = inner_case_outer_radius + magnet_coil_thickness
    magnet_coil_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=blanket_vessel_major_radius,b=magnet_coil_outer_radius*elongation*0.8,c=magnet_coil_outer_radius)
    outer_case_outer_radius = magnet_coil_outer_radius + outer_magnet_case_thickness
    outer_case_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=blanket_vessel_major_radius,b=outer_case_outer_radius*elongation*0.8,c=outer_case_outer_radius, boundary_type='vacuum')

    firstcm_outer_radius = inner_case_outer_radius+1
    firstcm_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=blanket_vessel_major_radius,b=firstcm_outer_radius*elongation*0.8,c=firstcm_outer_radius)
    
    midplane_upper_bound = openmc.ZPlane(z0=MIDPLANE_OFFSET)
    midplane_lower_bound = openmc.ZPlane(z0=-MIDPLANE_OFFSET)
    midplane_ib_bound = openmc.ZCylinder(r=blanket_vessel_major_radius)
    midpl = -midplane_upper_bound & +midplane_lower_bound & -midplane_ib_bound

    # Make two planes to cut the torus into a section
    # Angle follows right hand rule around z axis (https://www.desmos.com/3d/214a6bb908)
    section_angle_rad = np.radians(360*SECTION_CORRECTION)
    x_coeff, y_coeff = np.sin(section_angle_rad), -np.cos(section_angle_rad)
    xz_plane = openmc.Plane(a=0, b=1, boundary_type='periodic')
    angled_plane = openmc.Plane(a=x_coeff, b=y_coeff, boundary_type='periodic')
    xz_plane.periodic_surface = angled_plane
    torus_section = +xz_plane & +angled_plane

    # Create cells
    plasma_cell = openmc.Cell(
        name='plasma_cell',
        region=-plasma_surface & torus_section,
        fill=model_config['plasma_material']
    )

    first_wall_cell = openmc.Cell(
        name='first_wall_cell',
        region=+first_wall_inner_surface & -first_wall_outer_surface & torus_section,
        fill=model_config['first_wall_material']
    )
    first_wall_total_volume = total_layer_volume(R,first_wall_inner_radius,first_wall_outer_radius, elongation)
    if model_config['midplane_split']:
        first_wall_midpl = openmc.Cell(
            name='first_wall_midpl',
            region = first_wall_cell.region & midpl,
            fill=copy.deepcopy(model_config['first_wall_material'])
        )
        first_wall_cell.region = first_wall_cell.region & ~midpl
        first_wall_midpl_volume = peaking_sector_volume(R,first_wall_inner_radius,first_wall_outer_radius,first_wall_inner_radius*elongation,first_wall_outer_radius*elongation)
        first_wall_cell.fill.volume = first_wall_total_volume - first_wall_midpl_volume
        first_wall_midpl.fill.volume = first_wall_midpl_volume
    else:
        first_wall_cell.fill.volume = first_wall_total_volume

    cooling_channel_cell = openmc.Cell(
        name='cooling_channel_cell',
        region=+first_wall_outer_surface & -cooling_vessel_inner_surface & torus_section,
        fill=model_config['blanket_material']
    )
    cooling_channel_total_volume = total_layer_volume(R,first_wall_outer_radius,cooling_vessel_inner_radius, elongation)
    if model_config['midplane_split']:
        cooling_channel_midpl = openmc.Cell(
            name='cooling_channel_midpl',
            region = cooling_channel_cell.region & midpl,
            fill=copy.deepcopy(model_config['blanket_material'])
        )
        cooling_channel_cell.region = cooling_channel_cell.region & ~midpl
        cooling_channel_midpl_volume = peaking_sector_volume(R,first_wall_outer_radius,cooling_vessel_inner_radius,first_wall_outer_radius*elongation,cooling_vessel_inner_radius*elongation)
        cooling_channel_cell.fill.volume = cooling_channel_total_volume - cooling_channel_midpl_volume
        cooling_channel_midpl.fill.volume = cooling_channel_midpl_volume
    else:
        cooling_channel_cell.fill.volume = cooling_channel_total_volume

    cooling_vessel_cell = openmc.Cell(
        name='cooling_vessel_cell',
        region=+cooling_vessel_inner_surface & -cooling_vessel_outer_surface & torus_section,
        fill=model_config['cooling_vessel_material']
    )
    cooling_vessel_total_volume = total_layer_volume(R,cooling_vessel_inner_radius,cooling_vessel_outer_radius, elongation)
    if model_config['midplane_split']:
        cooling_vessel_midpl = openmc.Cell(
            name='cooling_vessel_midpl',
            region = cooling_vessel_cell.region & midpl,
            fill=copy.deepcopy(model_config['cooling_vessel_material'])
        )
        cooling_vessel_cell.region = cooling_vessel_cell.region & ~midpl
        cooling_vessel_midpl_volume = peaking_sector_volume(R,cooling_vessel_inner_radius,cooling_vessel_outer_radius,cooling_vessel_inner_radius*elongation,cooling_vessel_outer_radius*elongation)
        cooling_vessel_cell.fill.volume = cooling_vessel_total_volume - cooling_vessel_midpl_volume
        cooling_vessel_midpl.fill.volume = cooling_vessel_midpl_volume
    else:
        cooling_vessel_cell.fill.volume = cooling_vessel_total_volume

    vacuum_vessel_cell = openmc.Cell(
        name='vacuum_vessel_cell',
        region=+vacuum_vessel_inner_surface & -vacuum_vessel_outer_surface & torus_section,
        fill=model_config['vacuum_vessel_material']
    )
    vacuum_vessel_total_volume = total_layer_volume(R,vacuum_vessel_inner_radius,vacuum_vessel_outer_radius, elongation)
    if model_config['midplane_split']:
        vacuum_vessel_midpl = openmc.Cell(
            name='vacuum_vessel_midpl',
            region = vacuum_vessel_cell.region & midpl,
            fill=copy.deepcopy(model_config['vacuum_vessel_material'])
        )
        vacuum_vessel_cell.region = vacuum_vessel_cell.region & ~midpl
        vacuum_vessel_midpl_volume = peaking_sector_volume(R,vacuum_vessel_inner_radius,vacuum_vessel_outer_radius,vacuum_vessel_inner_radius*elongation,vacuum_vessel_outer_radius*elongation)
        vacuum_vessel_cell.fill.volume = vacuum_vessel_total_volume - vacuum_vessel_midpl_volume
        vacuum_vessel_midpl.fill.volume = vacuum_vessel_midpl_volume
    else:
        vacuum_vessel_cell.fill.volume = vacuum_vessel_total_volume

    # Volume of blanket is calculated differently because it has two non-concentric ellipses in its poloidal xs
    blanket_cell = openmc.Cell(
        name='blanket_cell',
        region=+vacuum_vessel_outer_surface & -blanket_vessel_inner_surface & torus_section,
        fill=model_config['blanket_material']
    )
    enclosed_volume = (2*np.pi*blanket_vessel_major_radius)*np.pi*blanket_vessel_inner_radius**2*elongation*0.8
    removed_volume = (2*np.pi*R)*np.pi*vacuum_vessel_outer_radius**2*elongation
    blanket_cell.fill.volume = (enclosed_volume - removed_volume) * SECTION_CORRECTION

    blanket_vessel_cell = openmc.Cell(
        name='blanket_vessel_cell',
        region=+blanket_vessel_inner_surface & -blanket_vessel_outer_surface & torus_section,
        fill=model_config['blanket_vessel_material']
    )
    total_blanket_vessel_volume = total_layer_volume(blanket_vessel_major_radius,blanket_vessel_inner_radius,blanket_vessel_outer_radius, elongation)
    if model_config['midplane_split']:
        blanket_vessel_midpl = openmc.Cell(
            name='blanket_vessel_midpl',
            region = blanket_vessel_cell.region & midpl,
            fill=copy.deepcopy(model_config['blanket_vessel_material'])
        )
        blanket_vessel_cell.region = blanket_vessel_cell.region & ~midpl
        blanket_vessel_midpl_volume = peaking_sector_volume(blanket_vessel_major_radius,blanket_vessel_inner_radius,blanket_vessel_outer_radius,blanket_vessel_inner_radius*elongation,blanket_vessel_outer_radius*elongation)
        blanket_vessel_cell.fill.volume = total_blanket_vessel_volume - blanket_vessel_midpl_volume
        blanket_vessel_midpl.fill.volume = blanket_vessel_midpl_volume
    else:
        blanket_vessel_cell.fill.volume = total_blanket_vessel_volume

    neutron_shield_cell = openmc.Cell(
        name='neutron_shield_cell',
        region=+blanket_vessel_outer_surface & -neutron_shield_outer_surface & torus_section,
        fill=model_config['magnetcase_material']
    )
    neutron_shield_cell.fill.volume = total_layer_volume(blanket_vessel_major_radius,blanket_vessel_outer_radius,neutron_shield_outer_radius, elongation)

    inner_case_cell = openmc.Cell(
        name='inner_case_cell',
        region=+neutron_shield_outer_surface & -inner_case_outer_surface & torus_section,
        fill=model_config['magnetcase_material']
    )
    inner_case_cell.fill.volume = total_layer_volume(blanket_vessel_major_radius,neutron_shield_outer_radius,inner_case_outer_radius, elongation)

    first_cm_cell = openmc.Cell(
        name='first_cm_cell',
        region=+inner_case_outer_surface & -firstcm_outer_surface & midpl & torus_section,
        fill=model_config['magnet_material']
    )
    first_cm_cell.fill.volume = peaking_sector_volume(blanket_vessel_major_radius,inner_case_outer_radius,firstcm_outer_radius,inner_case_outer_radius*elongation*0.8,firstcm_outer_radius*elongation*0.8)

    magnet_coil_cell = openmc.Cell(
        name='magnet_coil_cell',
        region=+inner_case_outer_surface & -magnet_coil_outer_surface & torus_section & ~first_cm_cell.region,
        fill=model_config['magnet_material']
    )
    magnet_coil_cell.fill.volume = total_layer_volume(blanket_vessel_major_radius,inner_case_outer_radius,magnet_coil_outer_radius, elongation) - \
        peaking_sector_volume(blanket_vessel_major_radius,inner_case_outer_radius,firstcm_outer_radius,inner_case_outer_radius*elongation*0.8,firstcm_outer_radius*elongation*0.8)

    outer_case_cell = openmc.Cell(
        name='outer_case_cell',
        region=+magnet_coil_outer_surface & -outer_case_outer_surface & torus_section,
        fill=model_config['magnetcase_material']
    )
    outer_case_cell.fill.volume = total_layer_volume(blanket_vessel_major_radius,magnet_coil_outer_radius,outer_case_outer_radius, elongation)

    universe = openmc.Universe(
        cells=[
            plasma_cell,
            first_wall_cell,
            cooling_channel_cell,
            cooling_vessel_cell,
            vacuum_vessel_cell,
            blanket_cell,
            blanket_vessel_cell,
            neutron_shield_cell,
            inner_case_cell,
            first_cm_cell,
            magnet_coil_cell,
            outer_case_cell
        ]
    )
    if model_config['midplane_split']:
        universe.add_cells([
            first_wall_midpl,
            cooling_channel_midpl,
            cooling_vessel_midpl,
            vacuum_vessel_midpl,
            blanket_vessel_midpl])

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
    # Make statepoints every 10 batches, ensuring the final batch is always included
    statepoint_set = set([i for i in range(10, model_config['batches']+1, 10)])
    statepoint_set.add(model_config['batches'])
    settings.statepoint = {'batches': list(statepoint_set)}
    settings.output = {'tallies': True}

    # # create weight windows for global variance reduction
    # ww_mesh = openmc.CylindricalMesh()
    # ww_mesh.r_grid = np.linspace((blanket_vessel_major_radius - outer_case_outer_radius), (blanket_vessel_major_radius + outer_case_outer_radius), 100)
    # ww_mesh.z_grid = np.linspace((0 - outer_case_outer_radius*elongation*0.8),(outer_case_outer_radius*elongation*0.8),100)
    # ww_mesh.phi_grid = [0, section_angle_rad]
    # ww_mesh_filter = openmc.MeshFilter(ww_mesh,filter_id=99)
    
    # wwg = openmc.WeightWindowGenerator(ww_mesh,np.logspace(-3, 8, 12))
    # wwg.max_realizations = 10
    # wwg.update_interval = 2
    # wwg.update_parameters = {'ratio' : 5.0,
    #                          'threshold': 0.5,
    #                          'value' : 'mean'}
    # settings.weight_window_generators = wwg

    #####################
    ## Define Tallies  ##
    #####################

    first_wall_cell_filter = openmc.CellFilter([first_wall_cell])
    vacuum_vessel_cell_filter = openmc.CellFilter([vacuum_vessel_cell])
    cooling_channel_cell_filter = openmc.CellFilter([cooling_channel_cell])
    cooling_vessel_cell_filter = openmc.CellFilter([cooling_vessel_cell])
    blanket_cell_filter = openmc.CellFilter([blanket_cell])
    blanket_vessel_cell_filter = openmc.CellFilter([blanket_vessel_cell])
    if model_config['midplane_split']:
        first_wall_midpl_filter = openmc.CellFilter([first_wall_midpl])
        cooling_channel_midpl_filter = openmc.CellFilter([cooling_channel_midpl])
        cooling_vessel_midpl_filter = openmc.CellFilter([cooling_vessel_midpl])
        vacuum_vessel_midpl_filter = openmc.CellFilter([vacuum_vessel_midpl])
        blanket_vessel_midpl_filter = openmc.CellFilter([blanket_vessel_midpl])
    
    magnetmat_filter = openmc.MaterialFilter([model_config['magnet_material']])
    first_cm_filter = openmc.CellFilter([first_cm_cell])

    # Tallies for flux and TBR in flibe mixture
    flibe_cell_filter = openmc.CellFilter([cooling_channel_cell, blanket_cell])
    energy_filter = openmc.EnergyFilter(np.logspace(0,7)) # 1eV to 100MeV

    flux_tally = openmc.Tally(tally_id=1, name="flux_blanket")
    flux_tally.filters = [flibe_cell_filter,energy_filter]
    flux_tally.scores = ["flux"]

    tbr_tally = openmc.Tally(tally_id=2, name='TBR')
    tbr_tally.filters = [flibe_cell_filter]
    tbr_tally.scores = ['(n,Xt)']

    # Tallies for neutron power deposition in each layer
    if model_config['photon_transport'] is True:
        heating_tally_type = ['heating']
    elif model_config['photon_transport'] is False:
        heating_tally_type = ['heating-local']

    first_wall_heating_tally = openmc.Tally(tally_id=3, name='neutron_heating_first_wall')
    first_wall_heating_tally.filters = [first_wall_cell_filter]
    first_wall_heating_tally.scores = heating_tally_type

    cooling_channel_heating_tally = openmc.Tally(tally_id=4, name='neutron_heating_cooling_channel')
    cooling_channel_heating_tally.filters = [cooling_channel_cell_filter]
    cooling_channel_heating_tally.scores = heating_tally_type

    cooling_vessel_heating_tally = openmc.Tally(tally_id=5, name='neutron_heating_cooling_vessel')
    cooling_vessel_heating_tally.filters = [cooling_vessel_cell_filter]
    cooling_vessel_heating_tally.scores = heating_tally_type

    vacuum_vessel_heating_tally = openmc.Tally(tally_id=6, name='neutron_heating_vacuum_vessel')
    vacuum_vessel_heating_tally.filters = [vacuum_vessel_cell_filter]
    vacuum_vessel_heating_tally.scores = heating_tally_type

    blanket_heating_tally = openmc.Tally(tally_id=7, name='neutron_heating_blanket')
    blanket_heating_tally.filters = [blanket_cell_filter]
    blanket_heating_tally.scores = heating_tally_type
    
    blanket_vessel_heating_tally = openmc.Tally(tally_id=8, name='neutron_heating_blanket_vessel')
    blanket_vessel_heating_tally.filters = [blanket_vessel_cell_filter]
    blanket_vessel_heating_tally.scores = heating_tally_type

    if model_config['midplane_split']:
        first_wall_midplane_tally = openmc.Tally(tally_id=23, name='neutron_heating_first_wall_midplane')
        first_wall_midplane_tally.filters = [first_wall_midpl_filter]
        first_wall_midplane_tally.scores = heating_tally_type

        cooling_channel_midpl_tally = openmc.Tally(tally_id=24, name='neutron_heating_cooling_channel_midplane')
        cooling_channel_midpl_tally.filters = [cooling_channel_midpl_filter]
        cooling_channel_midpl_tally.scores = heating_tally_type

        cooling_vessel_midpl_tally = openmc.Tally(tally_id=25, name='neutron_heating_cooling_vessel_midplane')
        cooling_vessel_midpl_tally.filters = [cooling_vessel_midpl_filter]
        cooling_vessel_midpl_tally.scores = heating_tally_type

        vacuum_vessel_midpl_tally = openmc.Tally(tally_id=26, name='neutron_heating_vacuum_vessel_midplane')
        vacuum_vessel_midpl_tally.filters = [vacuum_vessel_midpl_filter]
        vacuum_vessel_midpl_tally.scores = heating_tally_type

        blanket_vessel_midpl_tally = openmc.Tally(tally_id=28, name='neutron_heating_blanket_vessel_midplane')
        blanket_vessel_midpl_tally.filters = [blanket_vessel_midpl_filter]
        blanket_vessel_midpl_tally.scores = heating_tally_type

    high_e_filter = openmc.EnergyFilter([100000, 20000000])

    magnet_ff_tally = openmc.Tally(tally_id=9, name='fast_flux_magnet_coil')
    magnet_ff_tally.filters = [magnetmat_filter,high_e_filter]
    magnet_ff_tally.scores = ['flux']

    peak_magnet_ff_tally = openmc.Tally(tally_id=10, name='fast_flux_first_cm')
    peak_magnet_ff_tally.filters = [first_cm_filter,high_e_filter]
    peak_magnet_ff_tally.scores = ['flux']

    magnet_heating_tally = openmc.Tally(tally_id=11, name='fast_flux_magnet_coil')
    magnet_heating_tally.filters = [magnetmat_filter]
    magnet_heating_tally.scores = heating_tally_type

    peak_magnet_heating_tally = openmc.Tally(tally_id=12, name='fast_flux_first_cm')
    peak_magnet_heating_tally.filters = [first_cm_filter]
    peak_magnet_heating_tally.scores = heating_tally_type


    tallies = openmc.Tallies([
        flux_tally,
        tbr_tally,
        first_wall_heating_tally,
        cooling_channel_heating_tally,
        cooling_vessel_heating_tally,
        vacuum_vessel_heating_tally,
        blanket_heating_tally,
        blanket_vessel_heating_tally,
        magnet_ff_tally,
        peak_magnet_ff_tally,
        magnet_heating_tally,
        peak_magnet_heating_tally
    ])

    if model_config['midplane_split']:
        tallies.append(first_wall_midplane_tally)
        tallies.append(cooling_channel_midpl_tally)
        tallies.append(cooling_vessel_midpl_tally)
        tallies.append(vacuum_vessel_midpl_tally)
        tallies.append(blanket_vessel_midpl_tally)

    model = openmc.model.Model(
        geometry=geometry,
        settings=settings,
        tallies=tallies
    )

    return model



    