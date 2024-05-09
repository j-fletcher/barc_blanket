import xml.etree.ElementTree as ET
import openmc

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
        materials[material_id] = material

    return materials

# Example usage:
materials = import_materials('barc_blanket/materials/full_tank_inventory.xml')