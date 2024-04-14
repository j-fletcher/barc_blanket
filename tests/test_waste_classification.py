import openmc
import pytest

from barc_blanket.materials.waste_classification import check_class_c, sum_of_fractions, separate_tritium, make_activity_volume_density

class TestCheckClassC:

    def test_obviously_high_level_waste(self):
        """Test a material that is obviously high level waste, in this case pure Tc99"""
        material = openmc.Material(name='high_level_waste')
        material.add_nuclide('Tc99', 1.0)
        material.set_density('g/cm3', 11.5)
        result = check_class_c(material)
        # Ensure the result is 'False' as this is high level waste and print error message if not
        assert result == False, "Expected Tc99 to not meet the NRC requirements for class C low-level waste"

    def test_obviously_low_level_waste(self):
        """Test a material that is obviously low level waste, in this case pure O16"""
        material = openmc.Material(name='low_level_waste')
        material.add_nuclide('O16', 1.0)
        material.set_density('g/cm3', 0.2)
        result = check_class_c(material)
        # Ensure the result is 'True' as this is low level waste and print error message if not
        assert result == True, "Expected O16 to meet the NRC requirements for class C low-level waste"

class TestSumOfFractions:

    def test_nrc_example(self):
        """Ensure that the sum of fractions is calculated correctly by comparing to the NRC example
        in Paragraph 7 https://www.nrc.gov/reading-rm/doc-collections/cfr/part061/part061-0055.html"""

        # Create a material with 50 Ci/m3 of Sr-90 and 22 Ci/m3 of Cs-137
        nrc_example = {
            'Sr90': 50,
            'Cs137': 22
        }
        material = make_activity_volume_density(nrc_example)

        # Calculate the sum of fractions for column 1 of table 2
        sum_of_fractions_result, _ = sum_of_fractions(material, 2, 1)
        # Ensure the sum of fractions is greater than 1, indicating it is not class-A waste
        assert sum_of_fractions_result > 1, f"Expected sum of fractions to be greater than 1 but got {sum_of_fractions_result:0.2f}"

        # Calculate the sum of fractions for column 2 of table 2
        sum_of_fractions_result, _ = sum_of_fractions(material, 2, 2)
        # Ensure the sum of fractions is about 0.83, indicating it is class-B waste
        assert sum_of_fractions_result == pytest.approx(0.83, rel=0.01), f"Expected sum of fractions to be about 0.83 but got {sum_of_fractions_result:0.2f}"

class TestMakeActivityVolumeDensity:

    def test_nrc_example(self):
        """Attempt to create the same material as the NRC example in 
        Paragraph 7 https://www.nrc.gov/reading-rm/doc-collections/cfr/part061/part061-0055.html"""

        # Define target material in Ci/m3
        target_material = {
            'Sr90': 50,
            'Cs137': 22
        }

        # Create a material with 50 Ci/m3 of Sr-90 and 22 Ci/m3 of Cs-137
        material = make_activity_volume_density(target_material)

        activities_Bq_per_cm3 = material.get_activity(by_nuclide=True, units='Bq/cm3')

        # Ensure the material has the correct activity volume density
        for nuclide, target_activity in target_material.items():
            activity_Ci_per_m3 = (activities_Bq_per_cm3[nuclide] / 3.7e10) * 1e6

            assert activity_Ci_per_m3 == pytest.approx(target_activity, rel=0.01), f"Expected {nuclide} to have an activity of {target_activity:0.2f} Ci/m3 but got {activity_Ci_per_m3:0.2f} Ci/m3"


        