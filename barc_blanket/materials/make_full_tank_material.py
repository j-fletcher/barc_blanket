import numpy as np
import pandas as pd 
import openmc
from create_waste_material.py import *
'''
#################################################################
takes in pandas data frame as input argument
ex) 
df = pd.read_csv('Tanks_Slurry_Inventory - all_tank_data.csv')
total_waste_inventory = full_tank_inventory_material(df) --> openmc mixed material object of every tank except supernatant 241-B-201

This function calls the create_waste_material function already
##############################################################
'''


def full_tank_inventory_material(data):
	WasteSiteIDs = data['WasteSiteId'].dropna().unique()
	materials = []
	tank_phase_vol_dict = {}
	tanks_to_ignore = ['241-B-201_Supernatant']

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




