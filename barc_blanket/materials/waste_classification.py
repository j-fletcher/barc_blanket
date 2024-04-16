import openmc

CURIES_PER_BECQUEREL = 1/3.7e10 # NRC uses curies, OpenMC uses becquerels
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

def separate_nuclides(original_material:openmc.Material, nuclide_removal_efficiencies:dict):
    """Remove nuclides from a material with a given efficiency and return it as a new material
    with the remaining nuclides adjusted to maintain the same number of atoms,
    assuming the volume changes are linear with respect to mass density and independent of nuclide
    
    Parameters:
    -----------
    material: openmc.Material
        The material to remove nuclides from
    nuclide_removal_efficiencies: dict
        A dictionary of nuclides and their removal efficiencies in percent of total atoms
        key = nuclide, value = efficiency (float between 0 and 1)

    Returns:
    --------
    new_material: openmc.Material
        A new material with the same contents as the original but with some nuclides removed,
        and remaining concentrations adjusted accordingly
    """

    # The basic idea is that we want to remove some nuclides from the material while leaving the rest alone
    # We will accomplish this by first finding the mass density of all nuclides, then decreasing some of them
    # Then we will re-create the material with a mixture of the remaining nuclides using wt%
    # And finally we will adjust the density of the material to take into account what was removed

    # Get the nuclides in the material
    nuclides = original_material.get_nuclides()
    original_total_mass_density = original_material.get_mass_density()

    # Remove the nuclides with the given efficiencies
    # While doing this, keep track of each nuclide's mass density,
    # the total remaining mass density (should be lower than original),
    # and how much of the volume was removed (assuming 1 cm3 of the original material)
    remaining_mass_densities = {}
    remaining_total_mass_density = 0
    removed_volume = 0
    for nuclide in nuclides:
        original_mass_density = original_material.get_mass_density(nuclide)
        if nuclide in nuclide_removal_efficiencies.keys():
            efficiency = nuclide_removal_efficiencies[nuclide]
            # Ensure that efficiency is actually between 0 and 1
            if efficiency < 0 or efficiency > 1:
                raise ValueError(f"Removal efficiency must be between 0 and 1, but got {efficiency} for {nuclide}")
            new_mass_density = original_mass_density * (1 - efficiency)

            # Calculate how much volume the original nuclide took up
            volume_fraction = original_mass_density / original_total_mass_density
            # Calculate how much volume was removed
            removed_volume += volume_fraction * efficiency

            remaining_mass_density = new_mass_density
        else:
            remaining_mass_density = original_mass_density
        
        remaining_mass_densities[nuclide] = remaining_mass_density
        remaining_total_mass_density += remaining_mass_density

    # Create new material
    new_material = openmc.Material()

    # Scale each mass density to account for the change in volume
    # Assuming the volume change is linear with respect to mass density
    new_mass_densities = {}
    new_total_mass_density = 0
    for nuclide in nuclides:
        new_mass_density = remaining_mass_densities[nuclide] / (1 - removed_volume)
        new_mass_densities[nuclide] = new_mass_density
        new_total_mass_density += new_mass_density

    #new_total_mass_density = original_total_mass_density / (1 - removed_volume)

    # Add the remaining nuclides to the new material
    for nuclide, mass_density in new_mass_densities.items():
        if mass_density > 0:
            weight_percent = (mass_density / new_total_mass_density)*100
            new_material.add_nuclide(nuclide, weight_percent, 'wo')

    new_material.set_density('g/cm3', new_total_mass_density)
    
    return new_material

def remove_tritium(material:openmc.Material, efficiency):
    """Remove tritium from a material with a given efficiency
    
    Parameters:
    -----------
    material: openmc.Material
        The material to remove tritium from
    efficiency: float
        The efficiency of tritium removal (between 0 and 1)
    
    Returns:
    --------
    new_material: openmc.Material
        A new material with the same contents as the original but with tritium removed,
        and remaining concentrations adjusted accordingly
    """

    return separate_nuclides(material, {'H3': efficiency})

def remove_flibe(material:openmc.Material, efficiency):
    """Remove FLiBe from a material with a given efficiency
    
    Parameters:
    -----------
    material: openmc.Material
        The material to remove FLiBe from
    efficiency: float
        The efficiency of FLiBe removal (between 0 and 1)
    
    Returns:
    --------
    new_material: openmc.Material
        A new material with the same contents as the original but with FLiBe removed,
        and remaining concentrations adjusted accordingly
    """

    flibe_removal_dict = {
        'F19': efficiency,
        'Li6': efficiency,
        'Li7': efficiency,
        'Be9': efficiency,
        'Be10': efficiency
    }

    return separate_nuclides(material, flibe_removal_dict)


def vitrify_waste(material:openmc.Material, weight_percent_ratio):
    """Vitrify a material by adding a certain amount of borosilicate glass"""
    


def make_activity_volume_density(nuclide_activities_Ci_per_m3:dict):
    """Create a material with the given nuclides and activity concentrations in Ci/m3
    
    Parameters:
    -----------
    nuclide_activities_Ci_per_m3: dict
        A dictionary of nuclides and their activity concentrations in Ci/m3

    Returns:
    --------
    material: openmc.Material
        The new material
    """

    # Create the material
    material = openmc.Material()
    # For each nuclide, find how many kg of it are needed to achieve the given activity in 1 m3
    required_masses = {}
    total_required_mass = 0
    for nuclide, target_activity_Ci_per_m3 in nuclide_activities_Ci_per_m3.items():
        # Get the target activity in Bq/m3
        target_activity_Bq_per_m3 = target_activity_Ci_per_m3 / CURIES_PER_BECQUEREL
        decay_constant = openmc.data.decay_constant(nuclide)
        atoms_per_m3 = target_activity_Bq_per_m3 / decay_constant
        atomic_weight = openmc.data.atomic_mass(nuclide)
        kg_of_nuclide = (atoms_per_m3 * atomic_weight) * KG_PER_AMU
        total_required_mass += kg_of_nuclide
        required_masses[nuclide] = kg_of_nuclide

    # For each nuclide, add it to the material using wt%
    for nuclide, kg_of_nuclide in required_masses.items():
        weight_percent = (kg_of_nuclide / total_required_mass)*100
        material.add_nuclide(nuclide, weight_percent, 'wo')

    # Assuming the material is 1 m3
    material.set_density('kg/m3', total_required_mass)

    return material

