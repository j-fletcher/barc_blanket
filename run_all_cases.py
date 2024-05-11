# Run every case we are interested in.
# To be looked at tomorrow morning
import os

from barc_blanket.utilities import working_directory
from barc_blanket.models.barc_model_final import make_model, mo
from barc_blanket.materials.blanket_depletion import run_coupled_depletion
from barc_blanket.models.materials import flibe, lid, pbli, burner_mixture

CASES = {
    'pure_flibe': {'blanket_material': flibe(),
                   'name': "Pure FLiBe"},
    'pure_lid': {'blanket_material': lid(),
                 "name": "Pure LiD"},
    'pure_pbli': {'blanket_material': pbli(),
                  "name": "Pure PbLi"},
    'waste_01_flibe': {'blanket_material': burner_mixture(0.01, flibe=flibe()),
                       "name": "FLiBe 1% Full Tank Inventory"},
    'waste_01_lid': {'blanket_material': burner_mixture(0.01, flibe=lid()),
                     "name": "LiD 1% Full Tank Inventory"},
    'waste_01_pbli': {'blanket_material': burner_mixture(0.01, flibe=pbli()),
                      "name": "PbLi 1% Full Tank Inventory"},
    'waste_05_flibe': {'blanket_material': burner_mixture(0.05, flibe=flibe()),
                       "name": "FLiBe 5% Full Tank Inventory"},
    'waste_05_lid': {'blanket_material': burner_mixture(0.05, flibe=lid()),
                     "name": "LiD 5% Full Tank Inventory"},
    'waste_05_pbli': {'blanket_material': burner_mixture(0.05, flibe=pbli()),
                      "name": "PbLi 5% Full Tank Inventory"},
    'waste_10_flibe': {'blanket_material': burner_mixture(0.10, flibe=flibe()),
                       "name": "FLiBe 10% Full Tank Inventory"},
    'waste_10_lid': {'blanket_material': burner_mixture(0.10, flibe=lid()),
                     "name": "LiD 10% Full Tank Inventory"},
    'waste_10_pbli': {'blanket_material': burner_mixture(0.10, flibe=pbli()),
                      "name": "PbLi 10% Full Tank Inventory"},
}

BATCHES = 100
PARTICLES = 1e6
PHOTON_TRANSPORT = False

def main():
    for case, config in CASES.items():
        # create a working directory for each case
        os.makedirs(f"depletion_results/{case}", exist_ok=True)
        with working_directory(f"depletion_results/{case}"):
            model_config = {"batches": BATCHES,
                            "particles": PARTICLES,
                            "photon_transport": PHOTON_TRANSPORT,
                            "blanket_material": config['blanket_material']}
            
            model = make_model(model_config)
            model.export_to_model_xml()

            fusion_power = 2.2  # GW
            timesteps_years = [10] * 10 # 10 year timesteps for 100 years

            run_coupled_depletion(model, timesteps_years, fusion_power)

if __name__ == "__main__":
    main()