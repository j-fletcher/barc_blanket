import numpy as np
import pandas as pd 
import openmc
from create_waste_material.py import *
'''
#################################################################
takes in 2 inputs 
input 1: pandas data frame
input 2: integer between 0-3
			0 = full tank inventory
			1 = full tank inventory minus Pu/Th/U
			2 = A blend of the sludge phases of every tank (including stable and radionuclides), 
				plus only the radionuclides from the other phases; but with Cs and Sr removed from the final mixture 
			3 = A blend of the sludge phases of every tank (including stable and radionuclides) plus the radionuclides 
				from the other phases; but Cs, Sr, U/Th/Pu removed from the final mixture

			Defaults to 0
ex) 
df = pd.read_csv('Tanks_Slurry_Inventory - all_tank_data.csv')
total_waste_inventory = full_tank_inventory_material(df,0) --> openmc mixed material object of every tank except supernatant 241-B-201

This function calls the create_waste_material function already
##############################################################
'''


def full_tank_inventory_material(data,material_mix=0):
	radionuclide_list = ['Ru106','Cd113_m1','Sb125','Sn126','I129','Cs134','Cs137','Ba137_m1','C14','Sm151','Eu152','Eu154','Eu155',
							 'Ra226','Ac227','Ac228','Ra228','Th228','Th229','Th230','Pa231','Th232','U232','U233','U234','U235','U236',
							 'Np237','Pu238','U238','Pu239','Pu240','Am241','Pu241','Cm242','Pu242','Am243','Cm243','Cm244','H3','Ni59',
							 'Co60','Ni63','Se79','Sr90','Y90','Zr93','Nb93_m1','Nb94','Tc99']
	
	sludge_types = ['Sludge (Liquid & Solid)','Sludge Interstitial Liquid','Sludge Solid']

	tanks_to_ignore = ['241-B-201_Supernatant']


	if material_mix == 0:
		WasteSiteIDs = data['WasteSiteId'].dropna().unique()
		materials = []
		tank_phase_vol_dict = {}
		

		for tank_ID in WasteSiteIDs:
			for phase in data.loc[data['WasteSiteId'] == tank_ID]['WastePhase'].unique():
				if abs(data.loc[(data['WasteSiteId'] == tank_ID) & (data['WastePhase'] == phase),['Mass (kg)']].values[:,0].sum()) > 1e-12:
					tank_name = tank_ID+'_'+phase
					if tank_name not in tanks_to_ignore:
						waste_phase_volume_in_tank = data.loc[(data['WasteSiteId'] == tank_ID)& (data['WastePhase'] == phase)]['WastePhase Volume (L)'].unique().sum()
						mat = create_waste_material(tank_ID,phase,tank_name)
						materials.append(mat)
						tank_phase_vol_dict[tank_name] = waste_phase_volume_in_tank

		total_vol = sum(tank_phase_vol_dict.values())
		volume_fractions = []

		for volume in tank_phase_vol_dict:
			vol_frac = tank_phase_vol_dict[volume]/total_vol
			volume_fractions.append(vol_frac)

		waste_phase_volume_fractions = dict(zip(tank_phase_vol_dict.keys(),volume_fractions))

		total_tank_contents = openmc.Material.mix_materials(materials,volume_fractions,'vo')

		return total_tank_contents

	if material_mix == 1:
		WasteSiteIDs = data['WasteSiteId'].dropna().unique()
		materials = []
		tank_phase_vol_dict = {}

		for tank_ID in WasteSiteIDs:
			for phase in data.loc[data['WasteSiteId'] == tank_ID]['WastePhase'].unique():
				if abs(data.loc[(data['WasteSiteId'] == tank_ID) & (data['WastePhase'] == phase),['Mass (kg)']].values[:,0].sum()) > 1e-12:
					tank_name = tank_ID+'_'+phase
					if tank_name not in tanks_to_ignore:
						waste_phase_volume_in_tank = data.loc[(data['WasteSiteId'] == tank_ID)& (data['WastePhase'] == phase)]['WastePhase Volume (L)'].unique().sum()
						mat = create_waste_material(tank_ID,phase,tank_name)
						materials.append(mat)
						tank_phase_vol_dict[tank_name] = waste_phase_volume_in_tank

		total_vol = sum(tank_phase_vol_dict.values())
		volume_fractions = []

		for volume in tank_phase_vol_dict:
			vol_frac = tank_phase_vol_dict[volume]/total_vol
			volume_fractions.append(vol_frac)

		waste_phase_volume_fractions = dict(zip(tank_phase_vol_dict.keys(),volume_fractions))

		total_tank_contents = openmc.Material.mix_materials(materials,volume_fractions,'vo')
		total_tank_contents.remove_element('U')
		total_tank_contents.remove_element('Th')
		total_tank_contents.remove_element('Pu')

		return total_tank_contents

	if material_mix == 2:
		WasteSiteIDs = data['WasteSiteId'].dropna().unique()
		materials = []
		tank_phase_vol_dict = {}
		
		for tank_ID in WasteSiteIDs:
			for phase in data.loc[data['WasteSiteId'] == tank_ID]['WastePhase'].unique():
				if phase in sludge_types:
					if abs(data.loc[(data['WasteSiteId'] == tank_ID) & (data['WastePhase'] == phase),['Mass (kg)']].values[:,0].sum()) > 1e-12:
						tank_name = tank_ID+'_'+phase
						if tank_name not in tanks_to_ignore:
							waste_phase_volume_in_tank = data.loc[(data['WasteSiteId'] == tank_ID)& (data['WastePhase'] == phase)]['WastePhase Volume (L)'].unique().sum()
							mat = create_waste_material(tank_ID,phase,tank_name)
							materials.append(mat)
							tank_phase_vol_dict[tank_name] = waste_phase_volume_in_tank

				else:
					if abs(data.loc[(data['WasteSiteId'] == tank_ID) & (data['WastePhase'] == phase),['Mass (kg)']].values[:,0].sum()) > 1e-12:
						tank_name = tank_ID+'_'+phase
						if tank_name not in tanks_to_ignore:
							waste_phase_volume_in_tank = data.loc[(data['WasteSiteId'] == tank_ID)& (data['WastePhase'] == phase)]['WastePhase Volume (L)'].unique().sum()
							mat = create_waste_material(tank_ID,phase,tank_name)
							nuclides_in_mat = mat.get_nuclides()
							for nuclide in nuclides_in_mat:
								if nuclide not in radionuclide_list:
									mat.remove_nuclide(nuclide)
							materials.append(mat)
							tank_phase_vol_dict[tank_name] = waste_phase_volume_in_tank

		total_vol = sum(tank_phase_vol_dict.values())
		volume_fractions = []

		for volume in tank_phase_vol_dict:
			vol_frac = tank_phase_vol_dict[volume]/total_vol
			volume_fractions.append(vol_frac)

		waste_phase_volume_fractions = dict(zip(tank_phase_vol_dict.keys(),volume_fractions))

		total_tank_contents = openmc.Material.mix_materials(materials,volume_fractions,'vo')
		total_tank_contents.remove_element('Cs')
		total_tank_contents.remove_element('Sr')

		return total_tank_contents

	if material_mix == 3:
		WasteSiteIDs = data['WasteSiteId'].dropna().unique()
		materials = []
		tank_phase_vol_dict = {}
		
		for tank_ID in WasteSiteIDs:
			for phase in data.loc[data['WasteSiteId'] == tank_ID]['WastePhase'].unique():
				if phase in sludge_types:
					if abs(data.loc[(data['WasteSiteId'] == tank_ID) & (data['WastePhase'] == phase),['Mass (kg)']].values[:,0].sum()) > 1e-12:
						tank_name = tank_ID+'_'+phase
						if tank_name not in tanks_to_ignore:
							waste_phase_volume_in_tank = data.loc[(data['WasteSiteId'] == tank_ID)& (data['WastePhase'] == phase)]['WastePhase Volume (L)'].unique().sum()
							mat = create_waste_material(tank_ID,phase,tank_name)
							materials.append(mat)
							tank_phase_vol_dict[tank_name] = waste_phase_volume_in_tank

				else:
					if abs(data.loc[(data['WasteSiteId'] == tank_ID) & (data['WastePhase'] == phase),['Mass (kg)']].values[:,0].sum()) > 1e-12:
						tank_name = tank_ID+'_'+phase
						if tank_name not in tanks_to_ignore:
							waste_phase_volume_in_tank = data.loc[(data['WasteSiteId'] == tank_ID)& (data['WastePhase'] == phase)]['WastePhase Volume (L)'].unique().sum()
							mat = create_waste_material(tank_ID,phase,tank_name)
							nuclides_in_mat = mat.get_nuclides()
							for nuclide in nuclides_in_mat:
								if nuclide not in radionuclide_list:
									mat.remove_nuclide(nuclide)
							materials.append(mat)
							tank_phase_vol_dict[tank_name] = waste_phase_volume_in_tank

		total_vol = sum(tank_phase_vol_dict.values())
		volume_fractions = []

		for volume in tank_phase_vol_dict:
			vol_frac = tank_phase_vol_dict[volume]/total_vol
			volume_fractions.append(vol_frac)

		waste_phase_volume_fractions = dict(zip(tank_phase_vol_dict.keys(),volume_fractions))

		total_tank_contents = openmc.Material.mix_materials(materials,volume_fractions,'vo')
		total_tank_contents.remove_element('Cs')
		total_tank_contents.remove_element('Sr')
		total_tank_contents.remove_element('U')
		total_tank_contents.remove_element('Th')
		total_tank_contents.remove_element('Pu')

		return total_tank_contents








