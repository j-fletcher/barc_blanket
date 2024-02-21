# Plotting cross sections for arbitrary materials
# Author: Zander Keith [2024]

from openmc.material import Material
from openmc.plotter import plot_xs

def nu_fission(materials:list[Material], axis=None, normalize=True):
    """Plot the number of new neutrons produced per fission for a list of materials.
    
    Parameters:
    -----------
    materials: list[openmc.Material]
        A list of materials to plot the nu-fission for.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    normalize: bool
        If True, normalize the nu-fission cross section to the absorption cross section.
        If False, plot the nu-fission cross section as is.

    Returns:
    --------
    fig : matplotlib.figure.Figure
        If axis is None, then a Matplotlib Figure of the generated
        cross section will be returned. Otherwise, a value of
        None will be returned as the figure and axes have already been
        generated.

    """

    if normalize:
        divisor_xs = ['absorption']
    else:
        divisor_xs = None
    
    reactions = {}
    for material in materials:
        reactions[material] = ['nu-fission']

    fig = plot_xs(reactions, axis=axis, divisor_types=divisor_xs)

    return fig

def nu_scatter():
    pass

def absorption():
    pass