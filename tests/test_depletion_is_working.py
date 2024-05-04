import openmc
import openmc.deplete
import os
import numpy as np
import pytest
import matplotlib.pyplot as plt

from barc_blanket.utilities import working_directory

class TestCoupledDepletion:

    def test_Gd157_is_depleted(self):
        path = "tests/test_depletion_is_working/test_Gd157_is_depleted"
        # Clear the working directory
        os.system(f"rm -rf {path}")
        os.system(f"mkdir {path}")
        with working_directory(path):
            # create a simple model which is a torus section
            gd_material = openmc.Material(name="Gadolinium")
            gd_material.add_nuclide("Gd157", 1.0)
            gd_material.set_density("g/cm3", 7.9)
            gd_material.depletable = True

            R = 1.0
            inner_surface_radius = 0.2
            gadonlinium_width = 0.1
            outer_surface_radius = inner_surface_radius + gadonlinium_width

            inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=inner_surface_radius,c=inner_surface_radius)
            outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=outer_surface_radius,c=outer_surface_radius)

            bounding_sphere_surface = openmc.Sphere(r=R*2, boundary_type="vacuum")

            section_angle_rad = np.radians(45)
            x_coeff, y_coeff = np.sin(section_angle_rad), -np.cos(section_angle_rad)
            xz_plane = openmc.Plane(a=0, b=1, boundary_type='periodic')
            angled_plane = openmc.Plane(a=x_coeff, b=y_coeff, boundary_type='periodic')
            xz_plane.periodic_surface = angled_plane
            torus_section = +xz_plane & +angled_plane
            volume_correction = section_angle_rad/(2*np.pi)

            vacuum_cell = openmc.Cell(
                name="vacuum_cell",
                region=-inner_surface&torus_section,
                fill=None
            )
            
            gadolinium_cell = openmc.Cell(
                name="gadolinium_cell",
                region=+inner_surface&-outer_surface&torus_section,
                fill=gd_material
            )
            gadolinium_cell.fill.volume = np.pi * (outer_surface_radius**2 - inner_surface_radius**2) * volume_correction

            bounding_sphere_cell = openmc.Cell(
                name="bounding_sphere_cell",
                region=+outer_surface&-bounding_sphere_surface&torus_section,
                fill=None
            )

            universe = openmc.Universe()
            universe.add_cell(vacuum_cell)
            universe.add_cell(gadolinium_cell)
            universe.add_cell(bounding_sphere_cell)
            geometry = openmc.Geometry(universe)

            source = openmc.IndependentSource()
            source.particle = 'neutron'
            radius = openmc.stats.Discrete([R], [1])
            z_values = openmc.stats.Discrete([0], [1])
            angle = openmc.stats.Uniform(a=np.radians(0), b=section_angle_rad)
            source.space = openmc.stats.CylindricalIndependent(
                r=radius, phi=angle, z=z_values, origin=(0., 0., 0.))
            source.angle = openmc.stats.Isotropic() # Isotropic directio neutron is launched
            source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

            settings = openmc.Settings(run_mode='fixed source')
            settings.source=source
            settings.batches=10
            settings.inactive=2
            settings.particles=1000

            model=openmc.Model(geometry=geometry, 
                               settings=settings)

            timesteps_years = np.array([10, 10, 10, 10])
            timesteps_days = np.array(timesteps_years) * 365

            efus = 17.6e6  # eV
            ev2j = 1.60218e-19
            fusion_power = 2.2e9 # 2.2 GW
            neutron_rate = fusion_power / (efus * ev2j)  # n/s
            source_rates = np.ones_like(timesteps_years) * neutron_rate

            op = openmc.deplete.CoupledOperator(model, 
                                    reduce_chain=True, 
                                    reduce_chain_level=5, 
                                    normalization_mode='source-rate')
    
            openmc.deplete.CECMIntegrator(op, timesteps_days, source_rates=source_rates, timestep_units='d').integrate()

            # Load the depletion results
            results = openmc.deplete.Results(f"depletion_results.h5")

            material_at_times = []
            for i, time in enumerate(results.get_times()):
                materials = results.export_to_materials(burnup_index=i, path='materials.xml')
                material_at_times.append(materials[0])

            # Ensure the amount of Gd-157 is reduced at each time step
            atoms_at_times = [material.get_nuclide_atoms()['Gd157'] for material in material_at_times]

            for i in range(1, len(atoms_at_times)):
                assert atoms_at_times[i] < atoms_at_times[i-1], f"Expected Gd157 to be depleted but got {atoms_at_times[i]:0.2e}"

    def test_C14_is_depleted(self):
        path = "tests/test_depletion_is_working/test_C14_is_depleted"
        # Clear the working directory
        os.system(f"rm -rf {path}")
        os.system(f"mkdir {path}")
        with working_directory(path):
            # create a simple model which is a torus section
            c_material = openmc.Material(name="Carbon")
            c_material.add_nuclide("C14", 1.0)
            c_material.set_density("g/cm3", 2.2)
            c_material.depletable = True

            R = 1.0
            inner_surface_radius = 0.2
            gadonlinium_width = 0.1
            outer_surface_radius = inner_surface_radius + gadonlinium_width

            inner_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=inner_surface_radius,c=inner_surface_radius)
            outer_surface = openmc.ZTorus(x0=0,y0=0,z0=0,a=R,b=outer_surface_radius,c=outer_surface_radius)

            bounding_sphere_surface = openmc.Sphere(r=R*2, boundary_type="vacuum")

            section_angle_rad = np.radians(45)
            x_coeff, y_coeff = np.sin(section_angle_rad), -np.cos(section_angle_rad)
            xz_plane = openmc.Plane(a=0, b=1, boundary_type='periodic')
            angled_plane = openmc.Plane(a=x_coeff, b=y_coeff, boundary_type='periodic')
            xz_plane.periodic_surface = angled_plane
            torus_section = +xz_plane & +angled_plane
            volume_correction = section_angle_rad/(2*np.pi)

            vacuum_cell = openmc.Cell(
                name="vacuum_cell",
                region=-inner_surface&torus_section,
                fill=None
            )
            
            carbon_cell = openmc.Cell(
                name="carbon_cell",
                region=+inner_surface&-outer_surface&torus_section,
                fill=c_material
            )
            carbon_cell.fill.volume = np.pi * (outer_surface_radius**2 - inner_surface_radius**2) * volume_correction

            bounding_sphere_cell = openmc.Cell(
                name="bounding_sphere_cell",
                region=+outer_surface&-bounding_sphere_surface&torus_section,
                fill=None
            )

            universe = openmc.Universe()
            universe.add_cell(vacuum_cell)
            universe.add_cell(carbon_cell)
            universe.add_cell(bounding_sphere_cell)
            geometry = openmc.Geometry(universe)

            source = openmc.IndependentSource()
            source.particle = 'neutron'
            radius = openmc.stats.Discrete([R], [1])
            z_values = openmc.stats.Discrete([0], [1])
            angle = openmc.stats.Uniform(a=np.radians(0), b=section_angle_rad)
            source.space = openmc.stats.CylindricalIndependent(
                r=radius, phi=angle, z=z_values, origin=(0., 0., 0.))
            source.angle = openmc.stats.Isotropic() # Isotropic directio neutron is launched
            source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

            settings = openmc.Settings(run_mode='fixed source')
            settings.source=source
            settings.batches=10
            settings.inactive=2
            settings.particles=1000

            model=openmc.Model(geometry=geometry, 
                               settings=settings)

            timesteps_years = np.array([10, 10, 10, 10])
            timesteps_days = np.array(timesteps_years) * 365

            efus = 17.6e6  # eV
            ev2j = 1.60218e-19
            fusion_power = 2.2e9 # 2.2 GW
            neutron_rate = fusion_power / (efus * ev2j)  # n/s
            source_rates = np.ones_like(timesteps_years) * neutron_rate

            op = openmc.deplete.CoupledOperator(model, 
                                    reduce_chain=True, 
                                    reduce_chain_level=5, 
                                    normalization_mode='source-rate')
    
            openmc.deplete.CECMIntegrator(op, timesteps_days, source_rates=source_rates, timestep_units='d').integrate()

            # Load the depletion results
            results = openmc.deplete.Results(f"depletion_results.h5")

            material_at_times = []
            for i, time in enumerate(results.get_times()):
                materials = results.export_to_materials(burnup_index=i, path='materials.xml')
                material_at_times.append(materials[0])

            # Ensure the amount of C14 is reduced at each time step
            atoms_at_times = [material.get_nuclide_atoms()['C14'] for material in material_at_times]

            for i in range(1, len(atoms_at_times)):
                assert atoms_at_times[i] < atoms_at_times[i-1], f"Expected C14 to be depleted but got {atoms_at_times[i]:0.2e}"

    def test_multicell_material_depletion_just_works(self):
        """Verify that having one material in two different cells 'just works'

        Making a cylindrical setup divided into 4 parts.
        Upper right cell is A, lower right cell is B, lower left cell is C, upper left left cell is D.
        There is an isotropic neutron source in the center,
        and a shielding material in front of cells B and D.

        Cells A, B, C, and D will have the same material, but there will only be 3 material objects,
        so the material in cells C and D will be the same object.
        However, due to shielding each cell will have a distinct neutron flux and therefore different results.

        If having one material in multiple cells 'just works', then the rate of depletion
        for the material in cells C and D should be the average for the materials in A and B,
        because it is receiving a neutron flux that is a combination of the fluxes in A and B.
        """

        path = "tests/test_depletion_is_working/test_multicell_material_depletion_just_works"
        # Clear the working directory
        os.system(f"rm -rf {path}")
        os.system(f"mkdir {path}")

        with working_directory(path):

            batches = 20
            particles = 1000

            material_a = openmc.Material(name="Material A")
            material_b = openmc.Material(name="Material B")
            material_cd = openmc.Material(name="Material CD")

            materials = [material_a, material_b, material_cd]

            # Set the same nuclide in all materials
            for material in materials:
                material.add_nuclide("Gd157", 1.0)
                material.set_density("g/cm3", 1.0)
                material.depletable = True

            tungsten = openmc.Material(name='tungsten')
            tungsten.add_element('W', 1.0, 'wo')
            tungsten.set_density('g/cm3', 19.3)

            depletable_cell_outer_radius = 110
            depletable_cell_inner_radius = 100
            shield_cell_outer_radius = depletable_cell_inner_radius
            shield_cell_inner_radius = 90

            slab_height = 100

            # Create geometry
            top_plane = openmc.ZPlane(z0=slab_height/2, boundary_type='reflective')
            bottom_plane = openmc.ZPlane(z0=-slab_height/2, boundary_type='reflective')
            slab_region = -top_plane & +bottom_plane

            shield_surface = openmc.ZCylinder(r=shield_cell_inner_radius)
            depletable_inner_surface = openmc.ZCylinder(r=depletable_cell_inner_radius)
            depletable_outer_surface = openmc.ZCylinder(r=depletable_cell_outer_radius, boundary_type='vacuum')

            vacuum_region = -shield_surface & slab_region
            shield_region = +shield_surface & -depletable_inner_surface & slab_region
            depletable_region = +depletable_inner_surface & -depletable_outer_surface & slab_region

            x_plane = openmc.XPlane(x0=0)
            y_plane = openmc.YPlane(y0=0)

            vacuum_cell = openmc.Cell(name="Vacuum",
                                    region=vacuum_region,
                                    fill=None)

            cell_a_region = +x_plane & +y_plane & depletable_region
            shield_a_region = +x_plane & +y_plane & shield_region

            cell_b_region = -x_plane & +y_plane & depletable_region
            shield_b_region = -x_plane & +y_plane & shield_region

            cell_c_region = -x_plane & -y_plane & depletable_region
            shield_c_region = -x_plane & -y_plane & shield_region

            cell_d_region = +x_plane & -y_plane & depletable_region
            shield_d_region = +x_plane & -y_plane & shield_region

            cell_volume = np.pi * (depletable_cell_outer_radius**2 - depletable_cell_inner_radius**2) * slab_height / 4

            cell_a = openmc.Cell(name="Cell A",
                                region=cell_a_region,
                                fill=material_a)
            cell_b = openmc.Cell(name="Cell B",
                                region=cell_b_region,
                                fill=material_b)
            cell_c = openmc.Cell(name="Cell C",
                                region=cell_c_region,
                                fill=material_cd)
            cell_d = openmc.Cell(name="Cell D",
                                region=cell_d_region,
                                fill=material_cd)
            
            cell_a.fill.volume = cell_volume
            cell_b.fill.volume = cell_volume
            cell_c.fill.volume = 2*cell_volume

            shield_a = openmc.Cell(name="Shield A",
                                region=shield_a_region,
                                fill=None)
            shield_b = openmc.Cell(name="Shield B",
                                region=shield_b_region,
                                fill=tungsten)
            shield_c = openmc.Cell(name="Shield C",
                                region=shield_c_region,
                                fill=None)
            shield_d = openmc.Cell(name="Shield D",
                                region=shield_d_region,
                                fill=tungsten)

            universe = openmc.Universe(cells=[cell_a, cell_b, cell_c, cell_d,
                                            shield_a, shield_b, shield_c, shield_d,
                                            vacuum_cell])
            geometry = openmc.Geometry(universe)

            source = openmc.IndependentSource()
            source.particle = 'neutron'
            source.space = openmc.stats.Point((0,0,0))
            source.angle = openmc.stats.Isotropic() # Isotropic direction neutron is launched
            source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

            settings=openmc.Settings(run_mode='fixed source')
            settings.source = source
            settings.batches = batches
            settings.inactive = int(batches/5)
            settings.particles = particles

            model = openmc.model.Model(geometry=geometry, settings=settings)

            timesteps_days = np.array([10, 10, 10, 10])

            efus = 17.6e6  # eV
            ev2j = 1.60218e-19
            fusion_power = 2.2e9 # 2.2 GW
            neutron_rate = fusion_power / (efus * ev2j)  # n/s
            source_rates = np.ones_like(timesteps_days) * neutron_rate

            op = openmc.deplete.CoupledOperator(model, 
                                            reduce_chain=True, 
                                            reduce_chain_level=3, 
                                            normalization_mode='source-rate')
            
            openmc.deplete.CECMIntegrator(op, timesteps_days, source_rates=source_rates, timestep_units='d').integrate()

            # Load depletion results
            results = openmc.deplete.Results("depletion_results.h5")
            times = results.get_times()
            cell_a_results = results.get_atoms("1", "Gd157")[1]
            cell_b_results = results.get_atoms("2", "Gd157")[1]
            cell_cd_results = results.get_atoms("3", "Gd157")[1]

            # Plot the results and save to file (if you're interested)
            fig, ax = plt.subplots()
            ax.plot(times, cell_a_results, label="Cell A")
            ax.plot(times, cell_b_results, label="Cell B")
            ax.plot(times, cell_cd_results/2, label="Cell CD")
            ax.set_xlabel("Time [days]")
            ax.set_ylabel("Number of atoms")
            ax.legend()
            # Save fig as png
            plt.savefig("depletion_results.png")

            expected_final_cd = (cell_a_results[-1] + cell_b_results[-1])

            # Ensure the actual amount of Gd157 in cells C and D is close to the expected amount
            assert cell_cd_results[-1] == pytest.approx(expected_final_cd, rel=0.01), f"Expected Gd157 to have an activity of {expected_final_cd:0.2e} but got {cell_cd_results[-1]:0.2e}"

    def test_disconnected_cell_depletion_just_works(self):
        """Verify that having disconnected cells 'just works'

        Making a cylindrical setup divided into 4 parts.
        Upper right cel is A, lower left cell is C, and the remaining is cell BD
        There is an isotropic neutron source in the center,
        and a shielding material in front of cell A and half of BD

        Cells A, C, and BD will have the same material
        However, due to shielding only half of BD will have the same neutron flux as A and C,

        If having disconnected cells 'just works', then the rate of depletion
        for the material in cell BD should be the average for the materials in A and C,
        because it is receiving a neutron flux that is a combination of the fluxes in A and C.
        """

        path = "tests/test_depletion_is_working/test_disconnected_cell_depletion_just_works"
        # Clear the working directory
        os.system(f"rm -rf {path}")
        os.system(f"mkdir {path}")

        with working_directory(path):

            batches = 20
            particles = 1000

            material_a = openmc.Material(name="Material A")
            material_c = openmc.Material(name="Material C")
            material_bd = openmc.Material(name="Material BD")

            materials = [material_a, material_c, material_bd]

            # Set the same nuclide in all materials
            for material in materials:
                material.add_nuclide("Gd157", 1.0)
                material.set_density("g/cm3", 1.0)
                material.depletable = True

            tungsten = openmc.Material(name='tungsten')
            tungsten.add_element('W', 1.0, 'wo')
            tungsten.set_density('g/cm3', 19.3)

            depletable_cell_outer_radius = 110
            depletable_cell_inner_radius = 100
            shield_cell_inner_radius = 90

            slab_height = 100

            # Create geometry
            top_plane = openmc.ZPlane(z0=slab_height/2, boundary_type='reflective')
            bottom_plane = openmc.ZPlane(z0=-slab_height/2, boundary_type='reflective')
            slab_region = -top_plane & +bottom_plane

            shield_surface = openmc.ZCylinder(r=shield_cell_inner_radius)
            depletable_inner_surface = openmc.ZCylinder(r=depletable_cell_inner_radius)
            depletable_outer_surface = openmc.ZCylinder(r=depletable_cell_outer_radius, boundary_type='vacuum')

            vacuum_region = -shield_surface & slab_region
            shield_region = +shield_surface & -depletable_inner_surface & slab_region
            depletable_region = +depletable_inner_surface & -depletable_outer_surface & slab_region

            x_plane = openmc.XPlane(x0=0)
            y_plane = openmc.YPlane(y0=0)

            vacuum_cell = openmc.Cell(name="Vacuum",
                                    region=vacuum_region,
                                    fill=None)

            cell_a_region = +x_plane & +y_plane & depletable_region
            shield_a_region = +x_plane & +y_plane & shield_region

            cell_b_region = -x_plane & +y_plane & depletable_region
            shield_b_region = -x_plane & +y_plane & shield_region

            cell_c_region = -x_plane & -y_plane & depletable_region
            shield_c_region = -x_plane & -y_plane & shield_region

            cell_d_region = +x_plane & -y_plane & depletable_region
            shield_d_region = +x_plane & -y_plane & shield_region

            cell_volume = np.pi * (depletable_cell_outer_radius**2 - depletable_cell_inner_radius**2) * slab_height / 4

            cell_a = openmc.Cell(name="Cell A",
                                region=cell_a_region,
                                fill=material_a)
            cell_c = openmc.Cell(name="Cell C",
                                region=cell_c_region,
                                fill=material_c)
            cell_bd = openmc.Cell(name="Cell BD",
                                region=cell_b_region|cell_d_region,
                                fill=material_bd)
            
            cell_a.fill.volume = cell_volume
            cell_c.fill.volume = cell_volume
            cell_bd.fill.volume = 2*cell_volume

            shield_a = openmc.Cell(name="Shield A",
                                region=shield_a_region,
                                fill=tungsten)
            shield_b = openmc.Cell(name="Shield B",
                                region=shield_b_region,
                                fill=tungsten)
            shield_c = openmc.Cell(name="Shield C",
                                region=shield_c_region,
                                fill=None)
            shield_d = openmc.Cell(name="Shield D",
                                region=shield_d_region,
                                fill=None)

            universe = openmc.Universe(cells=[cell_a, cell_c, cell_bd,
                                            shield_a, shield_b, shield_c, shield_d,
                                            vacuum_cell])
            geometry = openmc.Geometry(universe)

            source = openmc.IndependentSource()
            source.particle = 'neutron'
            source.space = openmc.stats.Point((0,0,0))
            source.angle = openmc.stats.Isotropic() # Isotropic direction neutron is launched
            source.energy = openmc.stats.muir(e0=14.08e6, m_rat=5, kt=20000)

            settings=openmc.Settings(run_mode='fixed source')
            settings.source = source
            settings.batches = batches
            settings.inactive = int(batches/5)
            settings.particles = particles

            model = openmc.model.Model(geometry=geometry, settings=settings)

            timesteps_days = np.array([10, 10, 10, 10])

            efus = 17.6e6  # eV
            ev2j = 1.60218e-19
            fusion_power = 2.2e9 # 2.2 GW
            neutron_rate = fusion_power / (efus * ev2j)  # n/s
            source_rates = np.ones_like(timesteps_days) * neutron_rate

            op = openmc.deplete.CoupledOperator(model, 
                                            reduce_chain=True, 
                                            reduce_chain_level=3, 
                                            normalization_mode='source-rate')
            
            openmc.deplete.CECMIntegrator(op, timesteps_days, source_rates=source_rates, timestep_units='d').integrate()

            # Load depletion results
            results = openmc.deplete.Results("depletion_results.h5")
            times = results.get_times()
            cell_a_results = results.get_atoms("1", "Gd157")[1]
            cell_c_results = results.get_atoms("2", "Gd157")[1]
            cell_bd_results = results.get_atoms("3", "Gd157")[1]

            # Plot the results and save to file (if you're interested)
            fig, ax = plt.subplots()
            ax.plot(times, cell_a_results, label="Cell A")
            ax.plot(times, cell_c_results, label="Cell C")
            ax.plot(times, cell_bd_results/2, label="Cell BD")
            ax.set_xlabel("Time [days]")
            ax.set_ylabel("Number of atoms")
            ax.legend()
            # Save fig as png
            plt.savefig("depletion_results.png")

            expected_final_bd = (cell_a_results[-1] + cell_c_results[-1])

            # Ensure the actual amount of Gd157 in cells C and D is close to the expected amount
            assert cell_bd_results[-1] == pytest.approx(expected_final_bd, rel=0.01), f"Expected Gd157 to have an activity of {expected_final_bd:0.2e} but got {cell_bd_results[-1]:0.2e}"