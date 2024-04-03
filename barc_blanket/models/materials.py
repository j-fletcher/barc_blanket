import openmc

# Plasma
def dt_plasma():
    dt_plasma = openmc.Material(name='dt_plasma')
    dt_plasma.add_nuclide('H2', 1.0)
    dt_plasma.add_nuclide('H3', 1.0)
    dt_plasma.set_density('g/cm3', 1e-5)
    return dt_plasma

# FLIBE
def flibe(li6_enrichment=None):
    flibe = openmc.Material(name="flibe")
    flibe.depletable=True
    flibe.add_element("Be", 1.0, "ao")
    flibe.add_element("F", 4.0, "ao")
    flibe.set_density("g/cm3", 1.94)

    if li6_enrichment is None:
        flibe.add_element("Li", 2.0, "ao")
    else:
        flibe.add_element("Li", 2.0, "ao", 
                    enrichment=li6_enrichment, 
                    enrichment_target="Li6", 
                    enrichment_type="ao")

    return flibe

# Inconel 718 -
def inconel718():
    inconel718 = openmc.Material(name='inconel718')
    inconel718.depletable = True
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
    return inconel718

# Eurofer
def eurofer():
    eurofer = openmc.Material(name='eurofer')
    eurofer.depletable = True
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
    return eurofer

# V-4Cr-4Ti - pure -(from Segantin TRE https://github.com/SteSeg/tokamak_radiation_environment)
def v4cr4ti():
    v4cr4ti = openmc.Material(name='v4cr4ti')
    v4cr4ti.depletable = True
    v4cr4ti.add_element('V', 0.92, 'wo')
    v4cr4ti.add_element('Cr', 0.04, 'wo')
    v4cr4ti.add_element('Ti', 0.04, 'wo')
    v4cr4ti.set_density('g/cm3', 6.06)
    return v4cr4ti

# Tungsten - pure
def tungsten():
    tungsten = openmc.Material(name='tungsten')
    tungsten.depletable = True
    tungsten.add_element('W', 1.0, 'wo')
    tungsten.set_density('g/cm3', 19.3)
    return tungsten

# Water
def water():
    water = openmc.Material(name='water')
    water.depletable = True
    water.add_nuclide('H1', 2.0)
    water.add_nuclide('O16', 1.0)
    water.set_density('g/cm3', 1.0)
    return water

# Raw tank contents, do however you want to define this
# For now just say it's natural uranium TODO can change this to be whatever, or make a function for it
def simple_tank_contents():
    tank_contents = openmc.Material(name='tank_contents')
    tank_contents.depletable = True
    tank_contents.add_nuclide('U235', 0.0072)
    tank_contents.add_nuclide('U238', 0.9928)
    tank_contents.set_density('g/cm3', 19.1)
    return tank_contents

# Tank slurry contents
def burner_mixture(slurry_ratio, tank_contents=simple_tank_contents(), flibe=flibe()):
    """Create a mixture of flibe and tank contents for the burner blanket
    
    Parameters:
    ----------
    slurry_ratio : float
        The weight percent of slurry in the burner blanket
    tank_contents : openmc.Material, optional
        The tank contents to use in the mixture. Default is natural uranium.
    flibe : openmc.Material, optional
        The FLiBe material to use in the mixture. Default is the standard FLiBe material.
        Can pass in enriched flibe if desired

    Returns:
    -------
    burner_mixture : openmc.Material
        The mixture of FLiBe and tank contents
    
    """
    flibe_ao = 1 - slurry_ratio

    burner_mixture = openmc.Material.mix_materials(
        [flibe, tank_contents],
        [flibe_ao, slurry_ratio],
        'ao',
        name="burner_mixture"
    )
    burner_mixture.depletable = True

    return burner_mixture