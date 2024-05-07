from run_all_cases import CASES

from barc_blanket.utilities import working_directory
from barc_blanket.materials.blanket_depletion import postprocess_coupled_depletion, plot_results
from barc_blanket.models.barc_model_final import BLANKET_MATERIAL_ID

for case, config in CASES.items():
    try:
        with working_directory(f"depletion_results/{case}"):
            postprocess_coupled_depletion(BLANKET_MATERIAL_ID, remove_C14=False)
            plot_results(case, config['name'])
            postprocess_coupled_depletion(BLANKET_MATERIAL_ID, remove_C14=True)
            plot_results(f"{case}_no_C14", f"{config['name']} (no C14)")
    except:
        pass

