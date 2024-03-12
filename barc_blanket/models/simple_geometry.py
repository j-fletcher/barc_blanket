import openmc
import numpy as np

from .materials import dt_plasma, flibe, enriched_flibe, burner_mixture, v4cr4ti, tungsten

# Default model parameters
# TODO: this all assumes a circular cross-section, which is not necessarily the case
# Must determine if this is a reasonable assumption

DEFAULT_PARAMETERS = {
    'major_radius': 680,         # All dimensions are in cm
    'plasma_minor_radius': 120,
    'sol_width': 5,
    'first_wall_thickness': 1,   # How thick the plasma facing material is
    'vv_thickness': 2,           # How thick the wall of the vacuum vessel is
    'fusion_blanket_width': 20,  # Width of the material in the fusion blanket
    'burner_blanket_width': 100, # Width of the material in the burner blanket
    'bv_thickness': 5,           # How thick the burner vessel is

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

    #####################
    ## Assign Materials##
    #####################

    plasma_material = dt_plasma
    first_wall_material = tungsten
    vv_material = v4cr4ti
    fusion_blanket_material = enriched_flibe(model_config['li6_enrichment'])
    burner_blanket_material = burner_mixture(model_config['slurry_ratio'])
    bv_material = v4cr4ti

    
    #####################
    ## Define Geometry ##
    #####################

    R = model_config['major_radius']
    a = model_config['plasma_minor_radius']
    sol_width = model_config['sol_width']
    first_wall_thickness = model_config['first_wall_thickness']
    vv_thickness = model_config['vv_thickness']
    fusion_blanket_width = model_config['fusion_blanket_width']
    burner_blanket_width = model_config['burner_blanket_width']
    bv_thickness = model_config['bv_thickness']

    plasma_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=a,c=a)
    
    first_wall_inner_radius = a + sol_width
    first_wall_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=first_wall_inner_radius,c=first_wall_inner_radius)

    inboard_vv_inner_radius = first_wall_inner_radius + first_wall_thickness
    inboard_vv_outer_radius = inboard_vv_inner_radius + vv_thickness
    inboard_vv_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=inboard_vv_inner_radius,c=inboard_vv_inner_radius)
    inboard_vv_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=inboard_vv_outer_radius,c=inboard_vv_outer_radius)

    outboard_vv_inner_radius = inboard_vv_outer_radius + fusion_blanket_width
    outboard_vv_outer_radius = outboard_vv_inner_radius + vv_thickness
    outboard_vv_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=outboard_vv_inner_radius,c=outboard_vv_inner_radius)
    outboard_vv_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=outboard_vv_outer_radius,c=outboard_vv_outer_radius)

    bv_inner_radius = outboard_vv_outer_radius+burner_blanket_width
    bv_outer_radius = bv_inner_radius+bv_thickness
    bv_inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=bv_inner_radius,c=bv_inner_radius)
    bv_outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=bv_outer_radius,c=bv_outer_radius)

    bounding_sphere_surface = openmc.Sphere(r=2*R, boundary_type="vacuum")

    plasma_cell = openmc.Cell(
        name='plasma_cell',
        region=-plasma_surface,
        fill=plasma_material
    )

    sol_cell = openmc.Cell(
        name='sol_cell',
        region=+plasma_surface & -first_wall_inner_surface,
        fill=None
    )

    first_wall_cell = openmc.Cell(
        name='first_wall_cell',
        region=+first_wall_inner_surface & -inboard_vv_inner_surface,
        fill=first_wall_material
    )

    inboard_vv_cell = openmc.Cell(
        name='inboard_vv_cell',
        region=+inboard_vv_inner_surface & -inboard_vv_outer_surface,
        fill=vv_material
    )

    fusion_blanket_cell = openmc.Cell(
        name='fusion_blanket_cell',
        region=+inboard_vv_outer_surface & -outboard_vv_inner_surface,
        fill=fusion_blanket_material
    )

    outboard_vv_cell = openmc.Cell(
        name='outboard_vv_cell',
        region=+outboard_vv_inner_surface & -outboard_vv_outer_surface,
        fill=vv_material
    )

    burner_blanket_cell = openmc.Cell(
        name='burner_blanket_cell',
        region=+outboard_vv_outer_surface & -bv_inner_surface,
        fill=burner_blanket_material
    )

    bv_cell = openmc.Cell(
        name='bv_cell',
        region=+bv_inner_surface & -bv_outer_surface,
        fill=bv_material
    )

    bounding_sphere_cell = openmc.Cell(
        name='bounding_sphere_cell',
        region=+bv_outer_surface & -bounding_sphere_surface,
        fill=None
    )

    universe = openmc.Universe()
    universe.add_cell(plasma_cell)
    universe.add_cell(sol_cell)
    universe.add_cell(first_wall_cell)
    universe.add_cell(inboard_vv_cell)
    universe.add_cell(fusion_blanket_cell)
    universe.add_cell(outboard_vv_cell)
    universe.add_cell(burner_blanket_cell)
    universe.add_cell(bv_cell)
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