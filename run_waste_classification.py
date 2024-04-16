import openmc
import openmc.deplete
from barc_blanket.models.materials import burner_mixture
from barc_blanket.materials.waste_classification import check_class_c, sum_of_fractions, separate_nuclides, remove_flibe, remove_tritium, vitrify_waste

# Load tank contents from file

results = openmc.deplete.Results("./depletion_results_joe/depletion_results.h5")

times = results.get_times()

blanket_composition_at_time = []

for i, time in enumerate(times):
    materials = results.export_to_materials(burnup_index=i, path='./depletion_results_tank/materials.xml')
    blanket_composition_at_time.append(materials[5])

for original_material, time in zip(blanket_composition_at_time, times):

    #removed_tritium = remove_tritium(original_material, 0.9)
    #removed_flibe = remove_flibe(removed_tritium, 0.9)
    #sample_material = removed_flibe
    sample_material = original_material

    # Classify it as a waste material
    #result = check_class_c(sample_material)
    table_1_sum_of_fractions, _ = sum_of_fractions(sample_material, 1, None)
    table_2_sum_of_fractions, _ = sum_of_fractions(sample_material, 2, 3)

    print(f"Time: {time/365} years")
    print(f"Table 1 sum of fractions: {table_1_sum_of_fractions:0.2f}")
    print(f"Table 2 sum of fractions: {table_2_sum_of_fractions:0.2f}")
    #print(f"Is it class C waste? {result}")
