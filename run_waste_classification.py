from barc_blanket.models.materials import tank_contents
from barc_blanket.materials.waste_classification import check_class_c

# Make some tank contents

sample_material = tank_contents()

# Classify it as a waste material
result = check_class_c(sample_material)

print(result)
