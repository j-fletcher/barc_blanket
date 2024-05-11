import os
import openmc
import xml.etree.ElementTree as ET
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

def import_material_from_xml(materials_xml):
    """
    Create an openmc material by parsing an XML Materials file
    Maybe there's a better way of doing this
    """
    # Parse the materials.xml file
    tree = ET.parse(materials_xml)
    root = tree.getroot()
    print(root)
    # Create an empty dictionary to store materials
    materials = {}

    # Loop through each <material> tag in the XML file
    for material_elem in root.findall('material'):
        # Get material ID and name
        material_id = int(material_elem.get('id'))
        print(material_id)
        name = str(material_elem.get('name'))
        print(name)

        # Create a new OpenMC material object
        material = openmc.Material(material_id, name)

        # Loop through each <nuclide> tag within the material
        for nuclide_elem in material_elem.findall('nuclide'):
            nuclide = nuclide_elem.get('name')
            atom_percent = float(nuclide_elem.get('ao'))

            # Add the nuclide to the material with the specified atom percent
            material.add_nuclide(nuclide, atom_percent)

        # Store the material in the dictionary
        materials[str(material_id)] = material

    return materials

# Raw tank contents, do however you want to define this
def get_tank_contents(id=1):
    """
    ID information
    1 - most active tank, all waste phases present
    2 - All-tank mix, all waste phases present
    3 - All-tank mix, U, Th, Pu removed
    4 - All-tank mix, sludge phases only, Sr, Cs removed
    5 - All-tank mix, sludge phases only, Sr, Cs, U, Th, Pu removed
    6 - All-tank mix, sludge phases only
    7  - All-tank mix, sludge phases only, U, Th, Pu removed
    8 - All-tank mix, just sludge no radionuclides
    """
    if id == 1:
        # old workflow to get the most active tank material
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

    elif id == 2:
        material = import_material_from_xml("barc_blanket/materials/full_tank_inventory.xml")
        tank_contents = material["1371"]
    elif id == 3:
        material = import_material_from_xml("barc_blanket/materials/full_tank_inventory_no_PuThU.xml")
        tank_contents = material["1533"]
    elif id == 4:
        material = import_material_from_xml("barc_blanket/materials/sludge_plus_radionuclides_no_CsSr.xml")
        tank_contents = material["1059"]
    elif id == 5:
        material = import_material_from_xml("barc_blanket/materials/sludge_plus_radionuclides_no_CsSrPuThU.xml")
        tank_contents = material["1215"]
    elif id == 6:
        material = import_material_from_xml("barc_blanket/materials/full_tank_inventory_sludge_plus_radionuclides.xml")
        tank_contents = material["216"]
    elif id == 7:
        material = import_material_from_xml("barc_blanket/materials/full_tank_inventory_sludge_radionuclides_no_PuThU.xml")
        tank_contents = material["840"]
    elif id == 8:
        material = import_material_from_xml("barc_blanket/materials/full_tank_inventory_just_sludge.xml")
        tank_contents = material["60"]
    elif id == 9:
        tank_contents = openmc.Material()
        tank_contents.add_nuclide('Cs137', 0.5, 'wo')
        tank_contents.add_nuclide('Sr90', 0.5, 'wo')
        tank_contents.set_density('g/cm3', 2.3175)
    else:
        raise ValueError("Invalid ID for tank contents. See tank_contents function docstring for valid IDs.")

    tank_contents.depletable = True
    tank_contents.name = 'tank_contents'
    return tank_contents

# Mixture of tank contents and flibe for the blanket
def burner_mixture(slurry_ratio, id=1, flibe=flibe()):
    """Create a mixture of flibe and tank contents for the blanket
    
    Parameters:
    ----------
    slurry_ratio : float
        The weight percent of slurry in the blanket
    id : int
        Material ID to pass to the blanket mixture from tank_contents()
    flibe : openmc.Material, optional
        The FLiBe material to use in the mixture. Default is the standard FLiBe material.
        Can pass in enriched flibe if desired

    Returns:
    -------
    burner_mixture : openmc.Material
        The mixture of FLiBe and tank contents
    
    """
    tank_contents=get_tank_contents(id)

    flibe_ao = 1 - slurry_ratio

    burner_mixture = openmc.Material.mix_materials(
        [flibe, tank_contents],
        [flibe_ao, slurry_ratio],
        'wo',
        name="burner_mixture"
    )
    burner_mixture.depletable = True

    return burner_mixture