import openmc
from .create_waste import create_waste_material

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

    if li6_enrichment is None:
        flibe.add_element("Li", 2.0, "ao")
    else:
        flibe.add_element("Li", 2.0, "ao", 
                        enrichment=li6_enrichment, 
                        enrichment_target="Li6", 
                        enrichment_type="ao")

    flibe.set_density("g/cm3", 1.94)
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
def tank_contents():
    tank1_SCS = create_waste_material('241-A-106','Saltcake Solid','241-A-106-SCS')
    tank1_SLS = create_waste_material('241-A-106','Sludge (Liquid & Solid)','241-A-106-SLS')
    tank1_SS = create_waste_material('241-A-106','Sludge Solid','241-A-106-SS')

    #Tank Mixture waste phase volumes for 241-A-106 in L
    #SIL_vol = 198000
    SCS_vol = 84000
    SLS_vol = 38000+110000
    SS_vol = 41000

    total_vol = SCS_vol + SLS_vol + SS_vol

    #Waste Phase Volume Fraction
    #SIL_frac = SIL_vol/total_vol
    SCS_frac = SCS_vol/total_vol
    SLS_frac = SLS_vol/total_vol
    SS_frac = SS_vol/total_vol

    #Tank material object containing all present waste phases
    tank_contents = openmc.Material.mix_materials([tank1_SCS,tank1_SLS,tank1_SS],[SCS_frac,SLS_frac,SS_frac],'vo')

    tank_contents.depletable = True
    tank_contents.name = 'tank_contents'
    return tank_contents

# Mixture of tank contents and flibe for the blanket
def burner_mixture(slurry_ratio, tank_contents=tank_contents(), flibe=flibe()):
    """Create a mixture of flibe and tank contents for the blanket
    
    Parameters:
    ----------
    slurry_ratio : float
        The weight percent of slurry in the blanket
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
        'wo',
        name="burner_mixture"
    )
    burner_mixture.depletable = True

    return burner_mixture