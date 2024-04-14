import openmc
from barc_blanket.models.materials import tungsten, tank_contents
from barc_blanket.materials.waste_classification import check_class_c, sum_of_fractions

# Make some tank contents

#sample_material = tank_contents()
#sample_material = v4cr4ti()

# For quick test of tritium separation, looking at tritium in Technetium
# Tritium has no limits in class-c waste, but Technetium does
# So the idea is start with a very diluted mixture of Technetium in Tritium and make sure that 
# the class is NOT class C after removing all the Tritium

tritium = openmc.Material(name='tritium')
tritium.add_nuclide('H3', 1.0)
tritium.set_density('g/cm3', 0.1)

technetium = openmc.Material(name='technetium')
technetium.add_nuclide('Tc99', 1.0)
technetium.set_density('g/cm3', 11.5)

sample_material = openmc.Material.mix_materials([technetium, tritium], [0.001, 0.999], 'wo')

# Classify it as a waste material
result = check_class_c(sample_material)
table_1_sum_of_fractions, _ = sum_of_fractions(sample_material, 1, None)
table_2_sum_of_fractions, _ = sum_of_fractions(sample_material, 2, 3)

print(f"Table 1 sum of fractions: {table_1_sum_of_fractions:0.2f}")
print(f"Table 2 sum of fractions: {table_2_sum_of_fractions:0.2f}")
print(f"Is it class C waste? {result}")
