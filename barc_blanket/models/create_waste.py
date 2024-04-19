import numpy as np 
import pandas as pd 
import openmc
import warnings

def create_waste_material(tank,phase,mat_name):
    ############################### Compounds: Vectorize Data #######################################################################################
    C = openmc.data.atomic_weight('C')
    H = openmc.data.atomic_weight('H')
    O = openmc.data.atomic_weight('O')
    P = openmc.data.atomic_weight('P')
    N = openmc.data.atomic_weight('N')
    Cl = openmc.data.atomic_weight('Cl')
    F = openmc.data.atomic_weight('F')
    S = openmc.data.atomic_weight('S')
    elemmassvec = np.array([C,H,O,P,N,Cl,F,S])

    compounds = ["1-Butanol","1,1-Dichloroethene","1,1,1-Trichloroethane","1,1,2-Trichloro-1,2,2-trifluoroethane","1,1,2-Trichloroethane",
                 "1,1,2,2-Tetrachloroethane","1,2-Dichlorobenzene","1,2-Dichloroethane","1,2,4-Trichlorobenzene","1,4-Dichlorobenzene",
                 "2-Butanone","2-Chlorophenol","2-Ethoxyethanol","2-Methylphenol","2-Nitrophenol","2-Nitropropane","2,4-Dinitrotoluene",
                 "2,4,5-Trichlorophenol","2,4,6-Trichlorophenol","2,6-Bis(1,1-dimethylethyl)-4-methylphenol","4-Chloro-3-methylphenol",
                 "4-Methyl-2-Pentanone","4-Nitrophenol","Acenaphthene","Acetate","Acetone","Aroclors (Total PCB)","Benzene","Benzo(a)pyrene",
                 "Bis(2-ethylhexyl)phthalate","Butylbenzylphthalate","Carbon disulfide","Carbon tetrachloride","Chlorobenzene","Chloroform",
                 "CN","Cresol","Cresol (m & p)","Cyclohexanone","Di-n-butylphthalate","Di-n-octylphthalate","Dibenz[a,h]anthracene",
                 "Diethylphthalate","Diphenyl amine","Ethyl acetate","Ethyl ether","Ethylbenzene","Fluoranthene","Formate","Free OH","Glycolate",
                 "Hexachlorobenzene","Hexachlorobutadiene","Hexachloroethane","Hexone","Isobutanol","m-Cresol","Methylenechloride",
                 "Morpholine, 4-nitroso-","N-Nitroso-di-n-propylamine","N-Nitrosodimethylamine","Naphthalene","NH3","Nitrobenzene","NO2",
                 "NO3","Oxalate","Pentachlorophenol","Phenol","PO4","Pyrene","Pyridine","SO4","Sulfide","Tetrachloroethene","Thiosulfate",
                 "TIC as CO3","Toluene","Trans-1,3-Dichloropropene","Tributyl phosphate","Trichloroethene","Trichlorofluoromethane",
                 "Vinyl chloride","Xylene (m & p)","Xylene (o)","Xylenes (total)"]
    
    analytes_to_ignore = ['TOC','TotalAlpha','UTOTAL'] # these are accounted for in other surveys
    element_list = ['Ag','Al','As','B','Ba','Be','Bi','Br','Ca','Cd','Ce','Cl','Co','Cr','Cu','Eu','F','Fe','Hg','K','La','Li','Mg','Mn','Mo',
                    'Na','Nb','Nd','Ni','Pb','Pd','Pr','Rb','Rh','Ru','Sb','Se','Si','Sm','Sn','Sr','Ta','Te','Th','Ti','Tl','V','W','Y','Zn','Zr']
    radionuclide_list = ['106Ru','113mCd','125Sb','126Sn','129I','134Cs','137Cs','137mBa','14C','151Sm','152Eu','154Eu','155Eu','226Ra','227Ac',
                         '228Ac','228Ra','228Th','229Th','230Th','231Pa','232Th','232U','233U','234U','235U','236U','237Np','238Pu','238U',
                         '239Pu','240Pu','241Am','241Pu','242Cm','242Pu','243Am','243Cm','244Cm','3H','59Ni','60Co','63Ni','79Se','90Sr','90Y',
                         '93Zr','93mNb','94Nb','99Tc']

    formulae = ["C4H10O","C2H2Cl2","C2H3Cl3","C2Cl3F3","C2H3Cl3","C2H2Cl4","C6H4Cl2","C2H4Cl2","C6H3Cl3","C6H4Cl2","C4H8O","C6H5ClO","C4H10O2",
                "C7H8O","C6H5NO3","C3H7NO2","C7H6N2O4","C6H3Cl3O","C6H3Cl3O","C15H24O","C7H7ClO","C6H12O","C6H5NO3","C12H10","C2H3O2","C3H6O",
                "C12H3.5Cl6.5","C6H6","C20H12","C24H38O4","C19H20O4","CS2","CCl4","C6H5Cl","CHCl3","CN","C7H8O","C7H8O","C6H10O","C16H22O4",
                "C24H38O4","C22H14","C12H14O4","C12H11N","C4H8O2","C4H10O","C8H10","C16H10","CHO2","OH","C2H3O3","C6Cl6","C4Cl6","C2Cl6",
                "C6H12O","C4H10O","C7H8O","CH2Cl2","C4H8N2O2","C6H14N2O","C2H6N2O","C10H8","NH3","C6H5NO2","NO2","NO3","C2O4","C6HCl5O","C6H6O",
                "PO4","C16H10","C5H5N","SO4","S2","C2Cl4","O3S2","CO3","C7H8","C3H4Cl2","C12H27O4P","C2HCl3","CCl3F","C2H3Cl","C8H10","C8H10","C8H10"]
    
    ###################################[c,h,o,p,n,cl,f,s]
    vec_Butanol =                      [4,10,1,0,0,0,0,0]
    vec_Dichloroethene =               [2,2,0,0,0,2,0,0]
    vec_Trichloroethane1 =             [2,3,0,0,0,3,0,0]
    vec_Trichlorotrifluoroethane =     [2,0,0,0,0,3,3,0]
    vec_Trichloroethane2 =             [2,3,0,0,0,3,0,0]
    vec_Tetrachloroethane =            [2,2,0,0,0,4,0,0]
    vec_Dichlorobenzene2 =             [6,4,0,0,0,2,0,0]
    vec_Dichloroethane =               [2,4,0,0,0,2,0,0]
    vec_Trichlorobenzene =             [6,3,0,0,0,3,0,0]
    vec_Dichlorobenzene4 =             [6,4,0,0,0,2,0,0]
    vec_Butanone =                     [4,8,1,0,0,0,0,0]
    vec_Chlorophenol =                 [6,5,1,0,0,1,0,0]
    vec_Ethoxyethanol =                [4,10,2,0,0,0,0,0]
    vec_Methylphenol =                 [7,8,1,0,0,0,0,0]
    vec_Nitrophenol2 =                 [6,5,3,0,1,0,0,0]
    vec_Nitropropane =                 [3,7,2,0,1,0,0,0]
    vec_Dinitrotoluene =               [7,6,4,0,2,0,0,0]
    vec_Trichlorophenol5 =             [6,3,1,0,0,3,0,0]
    vec_Trichlorophenol6 =             [6,3,1,0,0,3,0,0]
    vec_Bisdimethylethylmethylphenol = [15,24,1,0,0,0,0,0]
    vec_Chloromethylphenol =           [7,7,1,0,0,1,0,0]
    vec_Methylpentanone =              [6,12,1,0,0,0,0,0]
    vec_Nitrophenol4 =                 [6,5,3,0,1,0,0,0]
    vec_Acenaphthene =                 [12,10,0,0,0,0,0,0]
    vec_Acetate =                      [2,3,2,0,0,0,0,0]
    vec_Acetone =                      [3,6,1,0,0,0,0,0]
    vec_Aroclors =                     [12,3.5,0,0,0,6.5,0,0]
    vec_Benzene =                      [6,6,0,0,0,0,0,0]
    vec_Benzopyrene =                  [20,12,0,0,0,0,0,0]
    vec_Bisethylhexylphthalate =       [24,38,4,0,0,0,0,0]
    vec_Butylbenzylphthalate =         [19,20,4,0,0,0,0,0]
    vec_Carbondisulfide =              [1,0,0,0,0,0,0,2]
    vec_Carbontetrachloride =          [1,0,0,0,0,4,0,0]
    vec_Chlorobenzene =                [6,5,0,0,0,1,0,0]
    vec_Chloroform =                   [1,1,0,0,0,3,0,0]
    vec_CN =                           [1,0,0,0,1,0,0,0]
    vec_Cresol =                       [7,8,1,0,0,0,0,0]
    vec_Cresolmp =                     [7,8,1,0,0,0,0,0]
    vec_Cyclohexanone =                [6,10,1,0,0,0,0,0]
    vec_Dibutylphthalate =             [16,22,4,0,0,0,0,0]
    vec_Dioctylphthalate =             [24,38,4,0,0,0,0,0]
    vec_Dibenzanthracene =             [22,14,0,0,0,0,0,0]
    vec_Diethylphthalate =             [12,14,4,0,0,0,0,0]
    vec_Diphenylamine =                [12,11,0,0,1,0,0,0]
    vec_Ethylacetate =                 [4,8,2,0,0,0,0,0]
    vec_Ethylether =                   [4,10,1,0,0,0,0,0]
    vec_Ethylbenzene =                 [8,10,0,0,0,0,0,0]
    vec_Fluoranthene =                 [16,10,0,0,0,0,0,0]
    vec_Formate =                      [1,1,2,0,0,0,0,0]
    vec_FreeOH =                       [0,1,1,0,0,0,0,0]
    vec_Glycolate =                    [2,3,3,0,0,0,0,0]
    vec_Hexachlorobenzene =            [6,0,0,0,0,6,0,0]
    vec_Hexachlorobutadiene =          [4,0,0,0,0,6,0,0]
    vec_Hexachloroethane =             [2,0,0,0,0,6,0,0]
    vec_Hexone =                       [6,12,1,0,0,0,0,0]
    vec_Isobutanol =                   [4,10,1,0,0,0,0,0]
    vec_mCresol =                      [7,8,1,0,0,0,0,0]
    vec_Methylenechloride =            [1,2,0,0,0,2,0,0]
    vec_Morpholinenitroso =            [4,8,2,0,2,0,0,0]
    vec_Nitrosodipropylamine =         [6,14,1,0,2,0,0,0]
    vec_Nitrosodimethylamine =         [0,0,0,0,2,0,0,0]
    vec_Naphthalene =                  [10,8,0,0,0,0,0,0]
    vec_NH3 =                          [0,3,0,0,1,0,0,0]
    vec_Nitrobenzene =                 [6,5,2,0,1,0,0,0]
    vec_NO2 =                          [0,0,2,0,1,0,0,0]
    vec_NO3 =                          [0,0,3,0,1,0,0,0]
    vec_Oxalate =                      [2,0,4,0,0,0,0,0]
    vec_Pentachlorophenol =            [6,1,1,0,0,5,0,0]
    vec_Phenol =                       [6,6,1,0,0,0,0,0]
    vec_PO4 =                          [0,0,4,1,0,0,0,0]
    vec_Pyrene =                       [16,10,0,0,0,0,0,0]
    vec_Pyridine =                     [5,5,0,0,1,0,0,0]
    vec_SO4 =                          [0,0,4,0,0,0,0,1]
    vec_Sulfide =                      [0,0,0,0,0,0,0,2]
    vec_Tetrachloroethene =            [2,0,0,0,0,4,0,0]
    vec_Thiosulfate =                  [0,0,3,0,0,0,0,2]
    vec_TICasCO3 =                     [1,0,3,0,0,0,0,0]
    vec_Toluene =                      [7,8,0,0,0,0,0,0]
    vec_Transdichloropropene =         [3,4,0,0,0,2,0,0]
    vec_Tributylphosphate =            [12,27,4,1,0,0,0,0]
    vec_Trichloroethene =              [2,1,0,0,0,3,0,0]
    vec_Trichlorofluoromethane =       [1,0,0,0,0,3,1,0]
    vec_Vinylchloride =                [2,3,0,0,0,1,0,0]
    vec_Xylenemp =                     [8,10,0,0,0,0,0,0]
    vec_Xyleneo =                      [8,10,0,0,0,0,0,0]
    vec_Xylenet =                      [8,10,0,0,0,0,0,0]

    mat_materials = np.stack((vec_Butanol,vec_Dichloroethene,vec_Trichloroethane1,vec_Trichlorotrifluoroethane,vec_Trichloroethane2,vec_Tetrachloroethane,
                              vec_Dichlorobenzene2,vec_Dichloroethane,vec_Trichlorobenzene,vec_Dichlorobenzene4,vec_Butanone,vec_Chlorophenol,
                              vec_Ethoxyethanol,vec_Methylphenol,vec_Nitrophenol2,vec_Nitropropane,vec_Dinitrotoluene,vec_Trichlorophenol5,
                              vec_Trichlorophenol6,vec_Bisdimethylethylmethylphenol,vec_Chloromethylphenol,vec_Methylpentanone,vec_Nitrophenol4,
                              vec_Acenaphthene,vec_Acetate,vec_Acetone,vec_Aroclors,vec_Benzene,vec_Benzopyrene,vec_Bisethylhexylphthalate,
                              vec_Butylbenzylphthalate,vec_Carbondisulfide,vec_Carbontetrachloride,vec_Chlorobenzene,vec_Chloroform,vec_CN,vec_Cresol,
                              vec_Cresolmp,vec_Cyclohexanone,vec_Dibutylphthalate,vec_Dioctylphthalate,vec_Dibenzanthracene,vec_Diethylphthalate,
                              vec_Diphenylamine,vec_Ethylacetate,vec_Ethylether,vec_Ethylbenzene,vec_Fluoranthene,vec_Formate,vec_FreeOH,vec_Glycolate,
                              vec_Hexachlorobenzene,vec_Hexachlorobutadiene,vec_Hexachloroethane,vec_Hexone,vec_Isobutanol,vec_mCresol,vec_Methylenechloride,
                              vec_Morpholinenitroso,vec_Nitrosodipropylamine,vec_Nitrosodimethylamine,vec_Naphthalene,vec_NH3,vec_Nitrobenzene,vec_NO2,
                              vec_NO3,vec_Oxalate,vec_Pentachlorophenol,vec_Phenol,vec_PO4,vec_Pyrene,vec_Pyridine,vec_SO4,vec_Sulfide,vec_Tetrachloroethene,
                              vec_Thiosulfate,vec_TICasCO3,vec_Toluene,vec_Transdichloropropene,vec_Tributylphosphate,vec_Trichloroethene,
                              vec_Trichlorofluoromethane,vec_Vinylchloride,vec_Xylenemp,vec_Xyleneo,vec_Xylenet))

    class Molecule:
        def __init__(self,name,formula,mass):
            self.name = name
            self.formula = formula
            self.mass = mass
        
        def add_elements(self,compvec):
            self._molecularmass = np.dot(compvec,elemmassvec)
            self.cfrac = compvec[0]*elemmassvec[0]/self._molecularmass
            self.hfrac = compvec[1]*elemmassvec[1]/self._molecularmass
            self.ofrac = compvec[2]*elemmassvec[2]/self._molecularmass
            self.pfrac = compvec[3]*elemmassvec[3]/self._molecularmass
            self.nfrac = compvec[4]*elemmassvec[4]/self._molecularmass
            self.clfrac = compvec[5]*elemmassvec[5]/self._molecularmass
            self.ffrac = compvec[6]*elemmassvec[6]/self._molecularmass
            self.sfrac = compvec[7]*elemmassvec[7]/self._molecularmass

    ############################### Compounds: Find and Decompose ###################################################################################
    # now correctly loops through compound library scan inside of analyte scan. if an analyte appears multiple times (i.e. in multiple waste types)
    # the mass will be summed across all instances of the analyte within the tank and phase of interest.
    # the summing is repeated and the number is overwritten for every instance of the analyte, but it never changes.
    def get_compound_masses_from_data(data,tankID,WastePhase):
        compound_masses = np.zeros(86)
        analytes = data.loc[(data['WasteSiteId'] == tankID) & (data['WastePhase'] == WastePhase),['Analyte']].values[:,0]
        for substance in analytes:
            if substance in compounds:
                compound_masses[compounds.index(substance)] = np.sum(data.loc[(data['WasteSiteId'] == tankID) & (data['WastePhase'] == WastePhase) & (data['Analyte'] == substance),['Mass (kg)']].values[:,0])
            elif substance == '144Ce/Pr' or substance == '239/240Pu' or substance == '243/244Cm':
                actvy = np.sum(data.loc[(data['WasteSiteId'] == tankID) & (data['WastePhase'] == WastePhase) & (data['Analyte'] == substance),['Activity (Ci)']].values[:,0])
                warnings.warn("Warning! Selected phase contains {} for which nuclide mass data cannot be determined! Activity present: {} Ci.".format(substance,actvy))
            elif substance not in element_list and substance not in radionuclide_list and substance not in analytes_to_ignore:
                raise KeyError('Unknown substance {} encountered in tank contents!'.format(substance))
    
        total_carbon_mass = 0
        total_hydrogen_mass = 0
        total_oxygen_mass = 0
        total_phosphorus_mass = 0
        total_nitrogen_mass = 0
        total_chlorine_mass = 0
        total_fluorine_mass = 0
        total_sulfur_mass = 0

        for com in compounds:
            compound = Molecule(com,formulae[compounds.index(com)],compound_masses[compounds.index(com)])
            compound.add_elements(mat_materials[compounds.index(com)])

            total_carbon_mass += compound.mass*compound.cfrac
            total_hydrogen_mass += compound.mass*compound.hfrac
            total_oxygen_mass += compound.mass*compound.ofrac
            total_phosphorus_mass += compound.mass*compound.pfrac
            total_nitrogen_mass += compound.mass*compound.nfrac
            total_chlorine_mass += compound.mass*compound.clfrac
            total_fluorine_mass += compound.mass*compound.ffrac
            total_sulfur_mass += compound.mass*compound.sfrac

        return [total_carbon_mass,total_hydrogen_mass,total_oxygen_mass,total_phosphorus_mass,
                total_nitrogen_mass,total_chlorine_mass,total_fluorine_mass,total_sulfur_mass]

    ############################### Elemental Surveys ###############################################################################################
    def element_masses_per_tank_per_waste_phase(data,tankID,WastePhase):
        # returns a dictionary of the surveyed elements and masses present
        # now correctly loops through element library scan inside of analyte scan. if an analyte appears multiple times (i.e. in multiple waste types)
        # the mass will be summed across all instances of the analyte within the tank and phase of interest
        analytes = df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase),['Analyte']].values[:,0]
        element_mass = []
        surveyed_elements = []

        for element in element_list:
            if element in analytes:
                element_mass.append(df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase)& (df['Analyte'] == element),['Mass (kg)']].values.sum())
                surveyed_elements.append(element)
        
        elements_dict = dict(zip(surveyed_elements,element_mass))

        return elements_dict

    ############################### Radionuclide Surveys ############################################################################################
    def nuclide_masses_per_tank_waste_phase(data,tankID,WastePhase):
        # returns a dictionary of the surveyed radionuclides and masses present
        # example: Cs137: 0.06
        # now correctly loops through nuclide library scan inside of analyte scan. if an analyte appears multiple times (i.e. in multiple waste types)
        # the mass will be summed across all instances of the analyte within the tank and phase of interest.
        analytes = df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase),['Analyte']].values[:,0]
        nuclides = []
        nuclides_reformatted = []
        masses = []
        
        for rn in radionuclide_list:
            if rn in analytes:
                nuclides.append(rn)
                masses.append(df.loc[(df['WasteSiteId'] == tankID) & (df['WastePhase'] == WastePhase)& (df['Analyte'] == rn),['Mass (kg)']].values.sum())
        nuclides_mass_dict = dict(zip(nuclides,masses))

        for nuclide in nuclides:
            digits = ''
            letters = ''
            new_key = ''
            for i in range(len(nuclide)):
                if nuclide[i].isdigit():
                    digits += nuclide[i]
                elif nuclide[i].isalpha():
                    if nuclide[i] == 'm':
                        if nuclide[i-1].isdigit():
                            digits += '_'
                            digits += nuclide[i]
                            digits += '1'
                        else:
                            letters += nuclide[i]
                    else:
                        letters += nuclide[i]
            new_key = letters+digits
            nuclides_reformatted.append(new_key)
        nuclides_reformatted_dict = dict(zip(nuclides_reformatted,nuclides_mass_dict.values()))
        return nuclides_reformatted_dict

    ############################### Remove Double-Counting Surveys ##################################################################################
    def elements_from_nuclides(data,tankID,WastePhase):
        # gives dictionary of elements which are represented in both radionuclide and total elemental surveys
        # form of output: dictionary of Element: Radionuclide
        # e.g. Cs: 137Cs
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

    def remove_duplicate_surveys(nuclides_present,elements_present,compounds_present):
        
        for element in ['C','H','O','P','N','Cl','F','S']:
            if element not in elements_present.keys():
                if compounds_present[element] > 0:
                    elements_present[element] = compounds_present[element]
        
        nuclide_double_counted_elements = []
        double_count_masses = {}

        for nuclide in nuclides_present.keys():
            letters = ''
            for i in range(len(nuclide)):
                if i < 2:
                    if nuclide[i].isalpha():
                        letters += nuclide[i]
            if letters in elements_present.keys():
                nuclide_double_counted_elements.append(letters)
                if letters in double_count_masses.keys():
                    double_count_masses[letters].append(nuclides_present[nuclide])
                else:
                    double_count_masses[letters] = []
                    double_count_masses[letters].append(nuclides_present[nuclide])

        nuclide_double_counted_elements = np.unique(nuclide_double_counted_elements)


        for element in nuclide_double_counted_elements:
            double_count_mass = np.sum(double_count_masses[element])
            elements_present[element] -= double_count_mass
        
        return [nuclides_present,elements_present]

    #################################################################################################################################################
    #                                                                                                                                               #
    #                               Parse the Data                                                                                                  #
    #                                                                                                                                               #
    #################################################################################################################################################
    df = pd.read_csv('/home/hallj/barc_blanket/barc_blanket/models/Tanks_Slurry_Inventory.csv')
    # check whether given phase is valid before proceeding
    phases = df.loc[(df['WasteSiteId'] == tank),['WastePhase']].values[:,0]
    if phase not in phases:
        raise KeyError('Waste phase {} not found in tank {}!'.format(phase,tank))
    
    # calculate masses of elements and radionuclides
    [cm,hm,om,pm,nm,clm,fm,sm] = get_compound_masses_from_data(df,tank,phase)
    present_compound_dict = dict(zip(['C','H','O','P','N','Cl','F','S'],[cm,hm,om,pm,nm,clm,fm,sm]))
    el_dict = element_masses_per_tank_per_waste_phase(df,tank,phase)
    rn_dict = nuclide_masses_per_tank_waste_phase(df,tank,phase)
    [all_radionuclides_present,all_elements_present] = remove_duplicate_surveys(rn_dict,el_dict,present_compound_dict)

    # Add Elements and Nuclides to Material (if mass above 1e-8 kg in parent tank & phase)
    # 1e-8 threshold is 11 mCi of Co-60, for example
    total_mass = np.sum(list(all_elements_present.values())) + np.sum(list(all_radionuclides_present.values()))
    waste_material = openmc.Material(name=mat_name)
    removals = []
    for el in all_elements_present.keys():
        if all_elements_present[el] > 1e-8:
            wf = all_elements_present[el] / total_mass
            waste_material.add_element(el,wf,'wo')
        else:
            removals.append(el)

    for rn in all_radionuclides_present.keys():
        if all_radionuclides_present[rn] > 1e-8:
            wf = all_radionuclides_present[rn] / total_mass
            waste_material.add_nuclide(rn,wf,'wo')
        else:
            removals.append(rn)
    if removals:
        warnings.warn("Warning! Removed the following elements/radionuclides from material as mass below 1e-8 kg:")
        print(removals)

    # Calculate density using a mass-weighted average across all waste types present in the phase of interest
    # this is still double-counting but as long as every analyte appears in every waste type it's fine 
    types = df.loc[(df['WasteSiteId'] == tank) & (df['WastePhase'] == phase),['WasteType']].values[:,0]
    utypes = np.unique(types)
    typemasses = np.zeros_like(utypes)
    typedensities = np.zeros_like(utypes)

    for t in utypes:
        tmass = df.loc[(df['WasteSiteId'] == tank) & (df['WastePhase'] == phase) & (df['WasteType'] == t),['Mass (kg)']].values[:,0]
        typemasses[np.where(utypes == t)] = np.nansum(tmass)
        # We ignore NaNs, which are found wherever the mixed nuclides (e.g. 239/240Pu) are listed
        tdens = df.loc[(df['WasteSiteId'] == tank) & (df['WastePhase'] == phase) & (df['WasteType'] == t),['ComponentDensity (g/mL)']].values[:,0]
        typedensities[np.where(utypes == t)] = np.nanmean(tdens)
    
    rho = np.average(typedensities,axis=None,weights=typemasses).round(3)
    waste_material.set_density('g/cm3',rho)

    return waste_material