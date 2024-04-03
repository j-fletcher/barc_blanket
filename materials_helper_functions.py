import numpy as np 
import pandas as pd 
import openmc

###############################################################################
'''
Takes in the total excel waste data, tank ID, desired waste phase, list of all compounds from total data, 
array of amount of individual molecules in each compound from total compounds list
outputs the mass of each element surveyed per tank per waste phase as a dictionary
'''
###############################################################################

def total_molecular_mass_per_molecule_dict(data,tankID,WastePhase,compounds,molecules):
    C = openmc.data.atomic_weight('C')
    H = openmc.data.atomic_weight('H')
    O = openmc.data.atomic_weight('O')
    P = openmc.data.atomic_weight('P')
    N = openmc.data.atomic_weight('N')
    Cl = openmc.data.atomic_weight('Cl')
    F = openmc.data.atomic_weight('F')
    S = openmc.data.atomic_weight('S')
    nuclides = ['C','H','O','P','N','Cl','F','S']
    elem_mass = np.array([C,H,O,P,N,Cl,F,S])
    compounds = compounds
    compound_mol_dict = dict(zip(compounds,molecules))
    total_compound_molecules = np.zeros((9))
    analytes = data.loc[(data['WasteSiteId'] == tankID) & (data['WastePhase'] == WastePhase),['Analyte'].values[:,0]
    for compound in analytes:
        if compound in compounds:
           total_compound_molecules += compound_mol_dict[compound]
    total_molecule_mass = total_compound_molecules*elem_mass
    nuclide_per_compound = dict(zip(nuclides,total_molecule_mass))
    return nuclide_per_compound

###############################################################################
'''
Takes in the total excel waste data, tank ID, desired waste phase, and the list of total individual elements surveyed and
outputs the mass of each element surveyed per tank per waste phase as a dictionary
'''
###############################################################################
def element_masses_per_tank_per_waste_phase(data,tankID,WastePhase,elements):
    analytes = df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase),['Analyte']].values[:,0]
    elements_dict = dict(zip(elements,np.zeros(len(elements))))
    for analyte in analytes:
        if analyte in elements:
            elements_dict[analyte] += df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase)
                                            & (df['Analyte'] == analyte),['Mass (kg)']].values[0,0]
    return elements_dict

##################################################################################
'''
Takes in the total excel waste data, tank ID, desired waste phase, and
outputs the mass of each indiviudally labled nuclide per tank per waste phase as a dictionary.
The nuclides are reformatted so that the element symbol is first followed by the mass number to make it easier to put into openmc.
'''
###################################################################################
def nuclide_masses_per_tank_waste_phase(data,tankID,WastePhase):
    analytes = df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase),['Analyte']].values[:,0]
    nuclides = []
    nuclides_reformatted = []
    masses = []
    for i in range(len(analytes)):
        if analytes[i][0].isdigit():
           if ((analytes[i][1].isdigit()) or ((analytes[i][1].isalpha()))):
               nuclides.append(analytes[i])
               masses.append(df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase)& (df['Analyte'] == analytes[i]),['Mass (kg)']].values[0,0])
    nuclides_mass_dict = dict(zip(nuclides,masses))
    for analyte in analytes:
        if analyte[0].isdigit():
            digits = ''
            letters = ''
            new_key = ''
            for i in range(len(analyte)):
                if analyte[0].isdigit():
                    if ((analyte[1].isdigit()) or ((analyte[1].isalpha()))):
                        if analyte[i].isdigit():
                                digits += analyte[i]
                        if analyte[i].isalpha():
                            if analyte[i] == 'm':
                                if analyte[i-1].isdigit():
                                    digits += analyte[i]
                                else:
                                    letters += analyte[i]
                            else:
                                letters += analyte[i]
            new_key = letters+digits
            nuclides_reformatted.append(new_key)
    nuclides_reformatted_dict = dict(zip(nuclides_reformatted,nuclides_mass_dict.values()))
    return nuclides_reformatted_dict
###########################################################################
'''
Get Compounds in Tank
'''
###########################################################################
def compounds_in_tank_list(data,tankID,WastePhase,compounds):
    analytes = df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase),['Analyte']].values[:,0]
    compounds_present = []
    for analyte in analytes:
        if analyte in compounds:
            compounds_present.append(analyte)
    return compounds_present

###########################################################################
'''
This Function takes creates a dictionary of the element name from a nuclide
ex. takes in {Co60 : 1.34e-9} and outputs {Co : Co60} 
'''
###########################################################################
def elements_from_nuclides(data,tankID,WastePhase):
    nuclides = nuclide_masses_per_tank_waste_phase(data,tankID,WastePhase)
    elements_separated = []
    for compounds in nuclides.keys():
        letters = ''
        for i in range(len(compounds)):
            sep_element = ''
            if i < 2:
                if compounds[i].isalpha():
                    letters += compounds[i]
        sep_element = letters
        elements_separated.append(sep_element)
        elements_separated_dict = dict(zip(elements_separated,nuclides.keys()))
    return elements_separated_dict
###############################################
'''
must get objects returned from:
    nuclide_masses_per_tank_waste_phase()
    element_masses_per_tank_per_waste_phase()
    total_molecular_mass_per_molecule_dict()
    elements_from_nuclides()
updates compound mass and element mass accordingly to avoid double counting
'''
################################################
def subtracting_nuclide_mass_from_comp_or_elements(nuclides_present,elements_present,compounds_present,elements_from_nuclides):
    double_elements = ['C','H','O','P','N','Cl','F','S']
    for substance in elements_present.keys():
        if substance in double_elements:
            if substance in elements_from_nuclides.keys():
                elements_present[substance] -= nuclides_present[elements_from_nuclides[substance]]
            elif compounds_present[substance] > 1e-8:
                compounds_present[substance] -= nuclides_present[elements_from_nuclides[substance]]
    return [elements_present,compounds_present]
