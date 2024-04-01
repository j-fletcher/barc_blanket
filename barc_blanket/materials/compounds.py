import openmc
import numpy as np
import math
import pandas as pd

################################################################## Constants: Atomic weights
c = 12.011
h = 1.008
o = 15.999
p = 30.973761998
n = 14.007
cl = 35.45
f = 18.998403162
s = 32.06

elemmassvec = np.array([c,h,o,p,n,cl,f,s])

################################################################## Vectorize data

compounds = ["1-Butanol","1,1-Dichloroethene","1,1,1-Trichloroethane","1,1,2-Trichloro-1,2,2-trifluoroethane","1,1,2-Trichloroethane",
             "1,1,2,2-Tetrachloroethane","1,2-Dichlorobenzene","1,2-Dichloroethane","1,2,4-Trichlorobenzene","1,4-Dichlorobenzene",
             "2-Butanone","2-Chlorophenol","2-Ethoxyethanol","2-Methylphenol","2-Nitrophenol","2-Nitropropane","2,4-Dinitrotoluene",
             "2,4,5-Trichlorophenol","2,4,6-Trichlorophenol","2,6-Bis(1,1-dimethylethyl)-4-methylphenol","4-Chloro-3-methylphenol",
             "4-Methyl-2-Pentanone","4-Nitrophenol","Acenaphthene","Acetate","Acetone","Aroclors(TotalPCB)","Benzene","Benzo(a)pyrene",
             "Bis(2-ethylhexyl)phthalate","Butylbenzylphthalate","Carbondisulfide","Carbontetrachloride","Chlorobenzene","Chloroform",
             "CN","Cresol","Cresol(m&p)","Cyclohexanone","Di-n-butylphthalate","Di-n-octylphthalate","Dibenz[a,h]anthracene",
             "Diethylphthalate","Diphenylamine","Ethylacetate","Ethylether","Ethylbenzene","Fluoranthene","Formate","FreeOH","Glycolate",
             "Hexachlorobenzene","Hexachlorobutadiene","Hexachloroethane","Hexone","Isobutanol","m-Cresol","Methylenechloride",
             "Morpholine,4-nitroso-","N-Nitroso-di-n-propylamine","N-Nitrosodimethylamine","Naphthalene","NH3","Nitrobenzene","NO2",
             "NO3","Oxalate","Pentachlorophenol","Phenol","PO4","Pyrene","Pyridine","SO4","Sulfide","Tetrachloroethene","Thiosulfate",
             "TICasCO3","Toluene","Trans-1,3-Dichloropropene","Tributylphosphate","Trichloroethene","Trichlorofluoromethane",
             "Vinylchloride","Xylene(m&p)","Xylene(o)","Xylenes(total)"]

formatted_comps = ["Butanol","Dichloroethene","Trichloroethane1","Trichlorotrifluoroethane","Trichloroethane2",
                    "Tetrachloroethane","Dichlorobenzene2","Dichloroethane","Trichlorobenzene","Dichlorobenzene4",
                    "Butanone","Chlorophenol","Ethoxyethanol","Methylphenol","Nitrophenol2","Nitropropane","Dinitrotoluene",
                    "Trichlorophenol5","Trichlorophenol6","Bisdimethylethylmethylphenol","Chloromethylphenol",
                    "Methylpentanone","Nitrophenol4","Acenaphthene","Acetate","Acetone","Aroclors","Benzene","Benzopyrene",
                    "Bisethylhexylphthalate","Butylbenzylphthalate","Carbondisulfide","Carbontetrachloride","Chlorobenzene","Chloroform",
                    "CN","Cresol","Cresolmp","Cyclohexanone","Dibutylphthalate","Dioctylphthalate","Dibenzanthracene",
                    "Diethylphthalate","Diphenylamine","Ethylacetate","Ethylether","Ethylbenzene","Fluoranthene","Formate","FreeOH","Glycolate",
                    "Hexachlorobenzene","Hexachlorobutadiene","Hexachloroethane","Hexone","Isobutanol","mCresol","Methylenechloride",
                    "Morpholinenitroso","Nitrosodipropylamine","Nitrosodimethylamine","Naphthalene","NH3","Nitrobenzene","NO2",
                    "NO3","Oxalate","Pentachlorophenol","Phenol","PO4","Pyrene","Pyridine","SO4","Sulfide","Tetrachloroethene","Thiosulfate",
                    "TICasCO3","Toluene","Transdichloropropene","Tributylphosphate","Trichloroethene","Trichlorofluoromethane",
                    "Vinylchloride","Xylenemp","Xyleneo","Xylenet"]

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

# example operation: total mass of each compound
# masses = mat_materials.dot(elemmassvec)

# associate each compound keyword with vectorized formula
compound_mol_dict = dict(zip(compounds,mat_materials))

################################################################## Molecule Operations

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
    
    def create_openmc_material(self):
        newmat = openmc.Material(name=self.name)
        newmat.add_elements_from_formula(self.formula)
        newmat.set_density('g/cm3',1) # dummy value
        self.openmc_material = newmat


################################################################## Read in Spreadsheet Data

def get_compound_masses_from_data(data,tankID,WastePhase):
    compound_masses = np.zeros(86)
    analytes = data.loc[(data['WasteSiteId'] == tankID) & (data['WastePhase'] == WastePhase),['Analyte']].values[:,0]
    for substance in compounds:
        if substance in analytes:
            compound_masses[compounds.index(substance)] = data.loc[(data['WasteSiteId'] == tankID) & (data['WastePhase'] == WastePhase) & (data['Analyte'] == substance),['Mass (kg)']].values[0,0]
        else:
            compound_masses[compounds.index(substance)] = 0
    
    ############################################################## Get Final Element Weight Fractions
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

################################################################## Verify mix_materials Method with Dummy Densities
# this is to check the method we'd discussed during the 2/23 meeting

#for com in compounds:
#    exec("%s = %d" % (formatted_comps[compounds.index(com)],
#         Molecule(com,formulae[compounds.index(com)],compound_wtfracs[compounds.index(com)])))
#    exec("%s.create_openmc_material()" % formatted_comps[compounds.index(com)])

## then mix materials into one, making sure to include the other isotopes that were added piecemeal
## get nuclide fractions
## compare to above method

################################################################# Example Calculation
df = pd.read_csv('Tanks_Slurry_Inventory - all_tank_data.csv')
[cm,hm,om,pm,nm,clm,fm,sm] = get_compound_masses_from_data(df,'241-C-103','Sludge (Liquid & Solid)')
print([cm,hm,om,pm,nm,clm,fm,sm])