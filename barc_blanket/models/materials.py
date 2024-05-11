import os
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

# Lithium deuteride
def lid():
    lid = openmc.Material(name='lid')
    lid.depletable = True
    lid.add_element('Li', 1.0, 'ao')
    lid.add_nuclide('H2', 1.0, 'ao')
    lid.set_density('g/cm3', 0.82)
    return lid

# Lead lithium
def pbli():
    pbli = openmc.Material(name='pbli')
    pbli.depletable = True
    pbli.add_element('Pb', 84.2, 'ao')
    pbli.add_element('Li', 15.8, 'ao')
    pbli.set_density('g/cm3', 10.2)
    return pbli

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

# SS316L for magnet case
def ss316L():
    ss316L = openmc.Material(name='ss316L')
    ss316L.add_element("Fe", 62.045, "wo")
    ss316L.add_element("Cr", 18., "wo")
    ss316L.add_element("Ni", 14., "wo")
    ss316L.add_element("Mo", 3., "wo")
    ss316L.add_element("Mn", 2., "wo")
    ss316L.add_element("Si", 0.75, "wo")
    ss316L.add_element("N", 0.1, "wo")
    ss316L.add_element("C", 0.03, "wo")
    ss316L.add_element("P", 0.045, "wo")
    ss316L.add_element("S", 0.03, "wo")
    ss316L.set_density("g/cm3", 7.99)

# Magnet winding pack mixture - combination of Stefano's material definition and Jack's
def magnetmat():
    copper = openmc.Material(name='copper')
    copper.add_element('Cu', 1.0)
    copper.set_density('g/cm3', 8.96)

    rebco = openmc.Material(name='rebco')
    rebco.add_element('Y',0.07734,'ao')
    rebco.add_element('Ba',0.154679,'ao')
    rebco.add_element('Cu',0.232019,'ao')
    rebco.add_element('O',0.535962,'ao')
    rebco.set_density('g/cm3',6.4)

    pbsn = openmc.Material(name='pbsn')
    pbsn.add_element('Pb', 0.37, percent_type='wo')
    pbsn.add_element('Sn', 0.63, percent_type='wo')
    pbsn.set_density('g/cm3', 8.8)

    hastelloy_c276 = openmc.Material(name='hastelloy_c276')
    hastelloy_c276.add_element('Ni', 0.5456, percent_type='wo')
    hastelloy_c276.add_element('Co', 0.025, percent_type='wo')
    hastelloy_c276.add_element('Cr', 0.16, percent_type='wo')
    hastelloy_c276.add_element('Mo', 0.16, percent_type='wo')
    hastelloy_c276.add_element('Fe', 0.05, percent_type='wo')
    hastelloy_c276.add_element('W', 0.04, percent_type='wo')
    hastelloy_c276.add_element('Mn', 0.01, percent_type='wo')
    hastelloy_c276.add_element('V', 0.0035, percent_type='wo')
    hastelloy_c276.add_element('Si', 0.0008, percent_type='wo')
    hastelloy_c276.add_element('C', 0.0001, percent_type='wo')
    hastelloy_c276.add_element('Cu', 0.005, percent_type='wo')
    hastelloy_c276.set_density('g/cm3', 8.89)

    silver = openmc.Material(name='silver')
    silver.add_element('Ag', 1.0)
    silver.set_density('g/cm3', 10.49)

    steel = ss316L()

    tape = openmc.Material.mix_materials([copper, silver, hastelloy_c276, rebco], [10/65.35, 3/65.35, 50/65.35, 2.35/65.35], percent_type='vo')
    magnetmat = openmc.Material.mix_materials([copper, steel, pbsn, tape], [0.15, 0.5, 0.05, 0.2], percent_type='vo', name="magnetmat")
    
    return magnetmat

# Raw tank contents, do however you want to define this
def tank_contents(mixture_name:str):
    """Return the material from the premade tank contents"""

    module_file_path = os.path.dirname(__file__)
    material_xml_path = f"{module_file_path}/../materials/{mixture_name}.xml"

    tank_contents = openmc.Materials.from_xml(material_xml_path)[0]

    return tank_contents

# Mixture of tank contents and flibe for the blanket
def burner_mixture(slurry_ratio, tank_contents=tank_contents("full_tank_inventory"), flibe=flibe()):
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