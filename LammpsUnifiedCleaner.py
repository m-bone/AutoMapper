import os
import re
import sys
from itertools import combinations_with_replacement
from LammpsTreatmentFuncs import clean_data, clean_settings, find_sections, get_data, add_section_keyword

# TO DO - Extras
# Remove lj/cut... for settings without hbonding. Would still need to fix hybrid in init

# ASSUMPTIONS
# Box size is identical between files

# That the header of the file operates in the standard format:
# Comment
#
#      a  atoms
#      b  bonds
#      c  angles
#      d  dihedrals
#      e  impropers

#      f  atom types
#      g bond types
#      h  angle types
#      i  dihedral types
#      j  improper types

#   x1 x2 xlo xhi
#   y1 y2 ylo yhi
#   z1 z2 zlo zhi

def file_unifier(directory, dataList, settingsFile):
    # Class for handling Lammps data
    class Data:
        def __init__(self, data):
            self.data = data
            self.sectionIndexList = find_sections(self.data)

            # Section data
            self.atoms = get_data('Atoms', self.data, self.sectionIndexList)
            self.masses = get_data('Masses', self.data, self.sectionIndexList)
            self.bonds = get_data('Bonds', self.data, self.sectionIndexList)
            self.angles = get_data('Angles', self.data, self.sectionIndexList)
            self.dihedrals = get_data('Dihedrals', self.data, self.sectionIndexList)
            self.impropers = get_data('Impropers', self.data, self.sectionIndexList)
        
        def get_atom_types(self):
            atom_types = {atom[2] for atom in self.atoms}
            return atom_types
        
        def get_bond_types(self):
            bond_types = {bond[1] for bond in self.bonds}
            return bond_types

        def get_angle_types(self):
            angle_types = {angle[1] for angle in self.angles}
            return angle_types

        def get_dihedral_types(self):
            dihedral_types = {dihedral[1] for dihedral in self.dihedrals}
            return dihedral_types

        def get_improper_types(self):
            improper_types = {improper[1] for improper in self.impropers}
            return improper_types

        def change_mass_types(self, unioned_atom_types):
            # Get masses that are using in atoms
            valid_masses = [mass for mass in self.masses if mass[0] in unioned_atom_types]
            # Get atom types
            mass_types = [mass[0] for mass in valid_masses]
            
            # Create dictionary of original atom type keys and new type values
            new_index_range = range(1, len(mass_types)+1)
            mass_zip = zip(mass_types, new_index_range)
            mass_change_dict = dict(mass_zip)
            
            # Change atom types to new types
            for massList in valid_masses:
                massList[0] = mass_change_dict[massList[0]]
            
            return valid_masses, mass_change_dict
        
        def change_atoms(self, mass_change_dict):
            for atomList in self.atoms:
                atomList[2] = mass_change_dict[atomList[2]]

    # Go to file directory
    os.chdir(directory)

    # Load files from dataList, tidy and initialise class object
    lammpsData = []
    for dataFile in dataList:
        with open(dataFile, 'r') as f:
            data = f.readlines()

        # Tidy data
        data = clean_data(data)

        # Initialise data class
        data = Data(data)
        lammpsData.append(data)


    # Load dataFile into python as a list of lists
    with open(settingsFile, 'r') as f:
        settings = f.readlines()
    
    # Tidy settings
    settings = clean_settings(settings)

    def union_types(typeAttr, lammpsData=lammpsData):
        lammpsTypes = []
        for data in lammpsData:
            # Get attribute for determining types
            func = getattr(data, typeAttr)
            # Run function and append to list
            lammpsTypes.append(func())

        # Union sets to remove duplicates and sort into numerical order list 
        types = sorted(set().union(*lammpsTypes))
        numTypes = len(types)
        
        print(f'{typeAttr}\n Types: {types}\n Count: {numTypes}')

        return types, numTypes

    # Union sets and create sorted list for each type
    atomTypes, numAtomTypes = union_types('get_atom_types')
    bondypes, numBondTypes = union_types('get_bond_types')
    angleTypes, numAngleTypes = union_types('get_angle_types')
    dihedralTypes, numdihedralTypes = union_types('get_dihedral_types')
    improperTypes, numImproperTypes = union_types('get_improper_types')

    # Build new masses section - eliminate unused masses
    processedMass, massDict = lammpsData[0].change_mass_types(atomTypes)
    processedMass = add_section_keyword('Masses', processedMass)

    # Update atom types and get atoms sections
    for data in lammpsData:
        data.change_atoms(massDict)
    processedAtoms = [add_section_keyword('Atoms', data.atoms) for data in lammpsData]

    # Update pair coeffs
    originalPairTuples = list(combinations_with_replacement(atomTypes, 2))

    print()

file_unifier('/home/matt/Documents/Bond_React_Python/', ['pre-system.data', 'post-system.data'], 'system.in.settings')