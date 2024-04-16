import openmc
from barc_blanket.models.materials import burner_mixture
from barc_blanket.materials.waste_classification import check_class_c, sum_of_fractions, separate_nuclides

# Make some tank contents

original_tank_contents = burner_mixture(slurry_ratio=0.05)

# Remove the Tritium and then FLiBe

TRITIUM_REMOVAL_EFFICIENCY = 0.95
FLIBE_REMOVAL_EFFICIENCY = 0.95
FLIBE_REMOVAL_DICT = {'F19': FLIBE_REMOVAL_EFFICIENCY, 
                      'Li6': FLIBE_REMOVAL_EFFICIENCY, 
                      'Li7': FLIBE_REMOVAL_EFFICIENCY, 
                      'Be9': FLIBE_REMOVAL_EFFICIENCY}

material_removed_tritium = separate_nuclides(original_tank_contents, {'H3': TRITIUM_REMOVAL_EFFICIENCY})
material_removed_flibe = separate_nuclides(material_removed_tritium, FLIBE_REMOVAL_DICT)

sample_material = material_removed_flibe

# Classify it as a waste material
result = check_class_c(sample_material)
table_1_sum_of_fractions, _ = sum_of_fractions(sample_material, 1, None)
table_2_sum_of_fractions, _ = sum_of_fractions(sample_material, 2, 3)

print(f"Table 1 sum of fractions: {table_1_sum_of_fractions:0.2f}")
print(f"Table 2 sum of fractions: {table_2_sum_of_fractions:0.2f}")
print(f"Is it class C waste? {result}")
