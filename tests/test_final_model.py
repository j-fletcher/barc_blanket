import os
import openmc.deplete

from barc_blanket.models.barc_model_final import make_model
from barc_blanket.utilities import working_directory
from barc_blanket.models.plot_geometry import plot_geometry

class TestFinalModel:

    def test_make_model(self):
        """Ensure the make_model function runs without error"""
        path = "tests/test_final_model/test_make_model"
        # Clear the working directory
        os.system(f"rm -rf {path}")
        os.system(f"mkdir {path}")

        with working_directory(path):

            try:
                model = make_model()
            except Exception as e:
                assert False, e

            # Plot the geometry of the model
            fig = plot_geometry(model)

            # Save the figure
            fig.savefig("geometry.png")

    def test_make_model_midplane_split(self):
        """Ensure the make_model function runs without error"""
        path = "tests/test_final_model/test_make_model_midplane_split"
        # Clear the working directory
        os.system(f"rm -rf {path}")
        os.system(f"mkdir {path}")

        with working_directory(path):

            try:
                model = make_model({'midplane_split': True})
            except Exception as e:
                assert False, e

            # Plot the geometry of the model
            fig = plot_geometry(model)

            # Save the figure
            fig.savefig("geometry.png")

            

