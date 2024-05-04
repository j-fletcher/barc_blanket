import openmc
import openmc.deplete
import pytest
import os
import numpy as np

from barc_blanket.utilities import working_directory, CROSS_SECTIONS, CHAIN_FILE

openmc.config['cross_sections'] = f"{CROSS_SECTIONS}"
openmc.config['chain_file'] = f"{CHAIN_FILE}"

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