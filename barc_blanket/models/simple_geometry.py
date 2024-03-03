# %%
import openmc
import numpy as np

# Default model parameters
# TODO: this all assumes a circular cross-section, which is not necessarily the case
# Must determine if this is a reasonable assumption

DEFAULT_PARAMETERS = {
    'major_radius': 680,         # All dimensions are in cm
    'plasma_minor_radius': 120,
    'sol_width': 5,
    'vv_thickness': 2,           # How thick the wall of the vacuum vessel is
    'fusion_blanket_width': 15,  # Width of the material in the fusion blanket
    'burner_blanket_width': 100, # Width of the material in the burner blanket
    'li6_enrichment': 0.076,     # atom% enrichment of Li6 in the FLiBe
    'slurry_ratio': 0.01         # wt% slurry in the burner blanket
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

    ######################
    ## Define Materials ##
    ######################

    # Plasma
    dt_plasma = openmc.Material(name='dt_plasma')
    dt_plasma.add_nuclide('H2', 1.0)
    dt_plasma.add_nuclide('H3', 1.0)
    dt_plasma.set_density('g/cm3', 1e-5)

    # FLIBE
    flibe = openmc.Material(name="flibe")
    flibe.add_element("Li", 2.0, "ao", 
                      enrichment=model_config['li6_enrichment'], 
                      enrichment_target="Li6", 
                      enrichment_type="ao")
    flibe.add_element("Be", 1.0, "ao")
    flibe.add_element("F", 4.0, "ao")
    flibe.set_density("g/cm3", 1.94)

    #TODO: add tank contents

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

    #####################
    ## Define Geometry ##
    #####################

    R = model_config['major_radius']
    a = model_config['plasma_minor_radius']
    sol_width = model_config['sol_width']
    vv_thickness = model_config['vv_thickness']
    fusion_blanket_width = model_config['fusion_blanket_width']
    burner_blanket_width = model_config['burner_blanket_width']

    plasma_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=a,c=a)

    inboard_vv_inner_radius = a + sol_width
    inboard_vv_outer_radius = inboard_vv_inner_radius + vv_thickness
    inboard_vv_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=inboard_vv_inner_radius,c=inboard_vv_inner_radius)
    inboard_vv_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=inboard_vv_outer_radius,c=inboard_vv_outer_radius)

    outboard_vv_inner_radius = inboard_vv_outer_radius + fusion_blanket_width
    outboard_vv_outer_radius = outboard_vv_inner_radius + vv_thickness
    outboard_vv_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=outboard_vv_inner_radius,c=outboard_vv_inner_radius)
    outboard_vv_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=outboard_vv_outer_radius,c=outboard_vv_outer_radius)

    burner_blanket_tank_inner_radius = outboard_vv_outer_radius+burner_blanket_width
    burner_blanket_tank_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=burner_blanket_tank_inner_radius,c=burner_blanket_tank_inner_radius)

    bounding_sphere_surface = openmc.Sphere(r=2*R, boundary_type="vacuum")

    plasma_cell = openmc.Cell(
        name='plasma_cell',
        region=-plasma_surface,
        fill=dt_plasma
    )

    sol_cell = openmc.Cell(
        name='sol_cell',
        region=+plasma_surface & -inboard_vv_inner_surface,
        fill=None
    )

    inboard_vv_cell = openmc.Cell(
        name='inboard_vv_cell',
        region=+inboard_vv_inner_surface & -inboard_vv_outer_surface,
        fill=v4cr4ti
    )

    fusion_blanket_cell = openmc.Cell(
        name='fusion_blanket_cell',
        region=+inboard_vv_outer_surface & -outboard_vv_inner_surface,
        fill=flibe
    )

    outboard_vv_cell = openmc.Cell(
        name='outboard_vv_cell',
        region=+outboard_vv_inner_surface & -outboard_vv_outer_surface,
        fill=v4cr4ti
    )

    # TODO: include slurry mixture
    burner_blanket_cell = openmc.Cell(
        name='burner_blanket_cell',
        region=+outboard_vv_outer_surface & -burner_blanket_tank_surface,
        fill=flibe
    )

    bounding_sphere_cell = openmc.Cell(
        name='bounding_sphere_cell',
        region=+burner_blanket_tank_surface & -bounding_sphere_surface,
        fill=None
    )

    universe = openmc.Universe()
    universe.add_cell(plasma_cell)
    universe.add_cell(sol_cell)
    universe.add_cell(inboard_vv_cell)
    universe.add_cell(fusion_blanket_cell)
    universe.add_cell(outboard_vv_cell)
    universe.add_cell(burner_blanket_cell)
    universe.add_cell(bounding_sphere_cell)
    geometry = openmc.Geometry(universe)

    #####################
    ## Define Settings ##
    #####################

    source = openmc.IndependentSource()
    source.particle = 'neutron'
    radius = openmc.stats.Discrete([R], [1]) # centered at major radius
    z_values = openmc.stats.Discrete([0], [1])
    angle = openmc.stats.Uniform(a=np.radians(0), b=np.radians(360))
    source.space = openmc.stats.CylindricalIndependent(
        r=radius, phi=angle, z=z_values, origin=(0., 0., 0.))
    source.angle = openmc.stats.Isotropic()
    source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

    settings = openmc.Settings(run_mode='fixed source')
    settings.photon_transport = False
    settings.source = source
    settings.batches = 50
    settings.particles = int(1e5) # modify this to shorten simulation, default was 1e6 
    settings.statepoint = {'batches': [
        5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100]}
    settings.output = {'tallies': True}

    #####################
    ## Define Tallies  ##
    #####################

    burner_cell_filter = openmc.CellFilter(burner_blanket_cell) # just burner blanket
    tbr_cell_filter = openmc.CellFilter([fusion_blanket_cell, burner_blanket_cell]) # includes both fusion and burner blankets
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

    tallies = openmc.Tallies([tally1,tally2,tally3]) 

    # Create model
    model = openmc.Model(geometry=geometry,
                         settings=settings, 
                         tallies=tallies)

    return model