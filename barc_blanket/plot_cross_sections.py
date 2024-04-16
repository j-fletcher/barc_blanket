# Plotting cross sections for arbitrary materials
# Author: Zander Keith [2024]

from numpy import inf, clip
from matplotlib import figure

from openmc.material import Material
from openmc.plotter import plot_xs

def clip_fig_y_axis(fig:figure.Figure, range=(1e-4, 1e4)):
    """ Sometimes the normalized cross sections can overflow, so clip  them to a reasonable range
    
    Parameters:
    -----------
    fig : matplotlib.figure.Figure
        The figure to clip the y-axis of.
    range : tuple
        The range to clip the y-axis to. Default is (1e-4, 1e4).
    """

    axes = fig.gca()
    min_y = inf
    max_y = -inf
    for line in axes.lines:
        clipped_data = clip(line.get_ydata(), range[0], range[1])
        min_y = min(min_y, min(clipped_data))
        max_y = max(max_y, max(clipped_data))

    axes.set_ylim(min_y, max_y)

##########################
# Neutron Multiplication #
##########################

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

    if normalize:
        clip_fig_y_axis(fig)

    return fig

def nu_scatter(materials:list[Material], axis=None, normalize=True):
    """Plot the number of new neutrons produced per scattering interaction for a list of materials.
    
    Parameters:
    -----------
    materials: list[openmc.Material]
        A list of materials to plot the nu-scatter for.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    normalize: bool
        If True, normalize the nu-scatter cross section to the absorption cross section.
        If False, plot the nu-scatter cross section as is.

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
        reactions[material] = ['nu-scatter']

    fig = plot_xs(reactions, axis=axis, divisor_types=divisor_xs)

    if normalize:
        clip_fig_y_axis(fig)

    return fig

def total_cross_section(materials:list[Material], axis=None, normalize=True):
    """Plot the total cross section for a list of materials.
    
    Parameters:
    -----------
    materials: list[openmc.Material]
        A list of materials to plot the total cross section for.
    axis : matplotlib.axes, optional
        A previously generated axis to use for plotting. If not specified,
        a new axis and figure will be generated.
    normalize: bool
        If True, normalize the total cross section to the absorption cross section.
        If False, plot the total cross section as is.

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
        reactions[material] = ['total']

    fig = plot_xs(reactions, axis=axis, divisor_types=divisor_xs)

    if normalize:
        clip_fig_y_axis(fig)

    return fig

def absorption():
    pass