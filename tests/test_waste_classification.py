import openmc

from barc_blanket.materials.waste_classification import check_class_c, sum_of_fractions, separate_tritium

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