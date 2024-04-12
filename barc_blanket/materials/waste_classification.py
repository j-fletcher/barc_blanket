import openmc

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

def check_class_c(material:openmc.Material, verbose=False):
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
    table_1_fractions: dict
        The relative fraction for each nuclide category in table 1
    table_2_fractions: dict
        The relative fraction for each nuclide category in table 2
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

    table_1_sum_of_fractions, table_1_fractions = sum_of_fractions(material, 1, None)
    if table_1_sum_of_fractions < 1 or verbose:
        table_2_sum_of_fractions, table_2_fractions = sum_of_fractions(material, 2, 3)
        if table_2_sum_of_fractions < 1:
            class_c = True
        else:
            class_c = False
    else:
        class_c = False

    if verbose:
        return class_c, table_1_fractions, table_2_fractions
    else:
        return class_c

