import openmc

CURIES_PER_BECQUEREL = 2.7e-11 # NRC uses curies, OpenMC uses becquerels
KG_PER_AMU = 1.66e-27
CUBIC_CENTIMETERS_PER_CUBIC_METER = 1e6

# Tables from https://www.nrc.gov/reading-rm/doc-collections/cfr/part061/part061-0055.html
# Assuming there is no 'activated metal' since it's a molten salt slurry

# Table 1: Concentration limits for waste classification in curies per Cubic Meter
TABLE_1_VOLUME_CONCENTRATION = {
    "C14": 8,
    "Tc99": 3,
    "I129": 0.08,
}

# Table 1 Concentration limits for waste classification in nanocuries per gram
TABLE_1_MASS_CONCENTRATION = {
    "long_lived_transuranic_alphas": 100,
    "Pu241": 3500,
    "Cm242": 20000,
}

# If something has 'no limit', put the value as None
TABLE_2_VOLUME_CONCENTRATION = {
    # Column 1
    1: {
        "all_short_lived_nuclides": 700,
        "H3": 40,
        "Co60": 700,
        "Ni63": 3.5,
        "Sr90": 0.04,
        "Cs137": 1
    },
    # Column 2
    2: {
        "all_short_lived_nuclides": None,
        "H3": None,
        "Co60": None,
        "Ni63": 70,
        "Sr90": 150,
        "Cs137": 44
    },
    # Column 3
    3: {
        "all_short_lived_nuclides": None,
        "H3": None,
        "Co60": None,
        "Ni63": 700,
        "Sr90": 7000,
        "Cs137": 4600
    }
}

def sum_of_fractions(material:openmc.Material, table, column):
    """Calculate the sum of fractions of a material
    See paragraph 7 on this page:
    https://www.nrc.gov/reading-rm/doc-collections/cfr/part061/part061-0055.html
    
    Parameters:
    -----------
    material: openmc.Material
        The material to check
    table: int
        The table to use for the calculation
    column: int
        The column to use for the calculation. Only valid for table 2.
    
    Returns:
    --------
    sum_of_fractions: float
        The sum of fractions for the material
    nuclide_fractions: dict
        The relative fraction for each nuclide category in the material
    """

    # Get the nuclides in the material
    nuclides = material.get_nuclides()

    # Determine which dictionary to use and fill in missing values if applicable
    if table == 1:
        volume_concentration = TABLE_1_VOLUME_CONCENTRATION
        mass_concentration = TABLE_1_MASS_CONCENTRATION

        for nuclide in nuclides:
            if nuclide not in mass_concentration.keys():
                # Check if it's an alpha-emitting transuranic with half life of greater than 5 years
                half_life_seconds = openmc.data.half_life(nuclide) 
                if half_life_seconds is not None: # If it's stable, this will be None
                    half_life_years = half_life_seconds / (365 * 24 * 60 * 60)
                    atomic_number = openmc.data.zam(nuclide)[0]
                    if atomic_number > 92 and half_life_years > 5:
                        # TODO: I'm pretty sure all unstable transuranic isotopes are alpha emitters,
                        # but we should double check this
                        mass_concentration[nuclide] = mass_concentration["long_lived_transuranic_alphas"]
    elif table == 2:
        if column is None:
            raise ValueError("Column must be specified for table 2")
        else:
            volume_concentration = TABLE_2_VOLUME_CONCENTRATION[column]
            mass_concentration = None

            for nuclide in nuclides:
                if nuclide not in volume_concentration.keys():
                    # Check if it has a half life of less than 5 years
                    half_life_seconds = openmc.data.half_life(nuclide)
                    if half_life_seconds is not None: # If it's stable, this will be None
                        half_life_years = half_life_seconds / (365 * 24 * 60 * 60)
                        if half_life_years < 5:
                            volume_concentration[nuclide] = volume_concentration["all_short_lived_nuclides"]
    else:
        raise ValueError("Invalid table number")

    # Get the activities in NRC units
    nuclide_activity_Bq_per_cm3 = material.get_activity(by_nuclide=True, units="Bq/cm3")
    nuclide_activity_Ci_per_m3 = {nuclide: activity * CURIES_PER_BECQUEREL * CUBIC_CENTIMETERS_PER_CUBIC_METER for nuclide, activity in nuclide_activity_Bq_per_cm3.items()}
    nuclide_activity_Bq_per_g = material.get_activity(by_nuclide=True, units="Bq/g")
    nuclide_activity_nCi_per_g = {nuclide: activity * CURIES_PER_BECQUEREL * 1e9 for nuclide, activity in nuclide_activity_Bq_per_g.items()}

    # Calculate the sum of fractions
    sum_of_fractions = 0
    nuclide_fractions = {}

    for nuclide in nuclides:
        if nuclide in volume_concentration.keys():
            if volume_concentration[nuclide] is not None:
                fraction = nuclide_activity_Ci_per_m3[nuclide] / volume_concentration[nuclide]
                sum_of_fractions += fraction
                nuclide_fractions[nuclide] = fraction
        elif mass_concentration is not None and nuclide in mass_concentration.keys():
            if mass_concentration[nuclide] is not None:
                fraction = nuclide_activity_nCi_per_g[nuclide] / mass_concentration[nuclide]
                sum_of_fractions += fraction
                nuclide_fractions[nuclide] = fraction

    return sum_of_fractions, nuclide_fractions

def check_class_c(material:openmc.Material):
    """Determine if the material is Class C waste according to the NRC

    https://www.nrc.gov/reading-rm/doc-collections/cfr/part061/part061-0055.html
    
    Parameters:
    -----------
    material: openmc.Material
        The material to check
    verbose: bool
        Return dictionary of results if True
    
    Returns:
    --------
    class_c: bool
        True if the material is Class C waste, False otherwise
    """

    # Stepping through the logic on the page linked above

    # Since we will likely be dealing with a mixture of long and short-lived waste,
    # Go to paragraph 5
    
    # (i) If the sum of fractions for table 1 does not exceed 0.1, 
    # the class shall be determined by the sum of fractions for table 2

    # (ii) If the sum of fractions for 1 exceeds 0.1, must ensure that it does not exceed
    # 1.0 for column 3 of table 2

    # Therefore, a material is class C low level waste if:
    # - The sum of fractions for table 1 does not exceed 1
    # - The sum of fractions for column 3 of table 2 does not exceed 1

    table_1_sum_of_fractions, _ = sum_of_fractions(material, 1, None)
    if table_1_sum_of_fractions < 1:
        table_2_sum_of_fractions, _ = sum_of_fractions(material, 2, 3)
        if table_2_sum_of_fractions < 1:
            class_c = True
        else:
            class_c = False
    else:
        class_c = False

    return class_c

def separate_tritium(original_material:openmc.Material, efficiency=1.0):
    """Remove tritium from a material with a given efficiency and return it as a new material
    
    Parameters:
    -----------
    material: openmc.Material
        The material to remove tritium from
    efficiency: float
        The efficiency of the separation process. Default is 1.0 (100%)

    Returns:
    --------
    new_material: openmc.Material
        A new material with the same contents as the original but with some tritium removed,
        and concentrations adjusted accordingly
    """

    # Get the nuclides in the material
    nuclides = original_material.get_nuclides()

    # Get the tritium concentration
    tritium_concentration = original_material.get_nuclide_atom_density("H3")

    # Calculate the new tritium concentration
    new_tritium_concentration = tritium_concentration * (1 - efficiency)

def make_activity_volume_density(nuclide_activity_Ci_per_m3:dict):
    """Create a material with the given nuclides and activity concentrations
    
    Parameters:
    -----------
    nuclide_activity_Ci_per_m3: dict
        A dictionary of nuclides and their activity concentrations in Ci/m3

    Returns:
    --------
    material: openmc.Material
        The new material
    """

    nuclides = nuclide_activity_Ci_per_m3.keys()

    # Create the material
    material = openmc.Material()
    # For each nuclide, add it to the material with weight% assuming there is 1kg of the material
    for nuclide in nuclides:
        # Get the target activity in Bq/m3
        target_activity_Bq_per_m3 = nuclide_activity_Ci_per_m3[nuclide] / CURIES_PER_BECQUEREL
        decay_constant = openmc.data.decay_constant(nuclide)
        atoms_per_m3 = target_activity_Bq_per_m3 / decay_constant
        atomic_weight = openmc.data.atomic_mass(nuclide)
        kg_of_nuclide = (atoms_per_m3 / atomic_weight) * KG_PER_AMU
        material.add_nuclide(nuclide, kg_of_nuclide, 'wo')

    material.set_density('kg/m3', 1)
    return material

