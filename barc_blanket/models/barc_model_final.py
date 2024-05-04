import openmc

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

    'section_angle': 360/14, # 14 coils -> this many degrees per coil

    'batches': 50,
    'inactive_batches': 5,
    'particles': int(1e5),

    'photon_transport': False
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

    R = model_config['major_radius']
    a = model_config['minor_radius']
    elongation = model_config['elongation']
    b=elongation*a

    first_wall_thickness = model_config['first_wall_thickness']
    cooling_channel_width = model_config['cooling_channel_width']
    cooling_vessel_thickness = model_config['cooling_vessel_thickness']
    vacuum_vessel_thickness = model_config['vacuum_vessel_thickness']
    blanket_vessel_thickness = model_config['blanket_vessel_thickness']

    # Taking the provided minor radius as the smaller one
    # This is gonna make the thicknesses on the top and bottom a little bigger,
    # but I don't think that will make a significant difference
    plasma_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=b,c=a)

    first_wall_inner_radius = a
    first_wall_outer_radius = a + first_wall_thickness
    first_wall_inner_surface = openmc.ZTorus()

    cooling_vessel_inner_radius = first_wall_outer_radius # These two are in contact
    cooling_vessel_outer_radius = cooling_vessel_inner_radius + cooling_vessel_thickness

    vacuum_vessel_inner_radius = cooling_vessel_outer_radius
    vacuum_vessel_outer_radius = vacuum_vessel_inner_radius + vacuum_vessel_thickness

    blanket_vessel_inner_radius = vacuum_vessel_outer_radius