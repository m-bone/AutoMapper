import os
import re
import sys
from itertools import combinations_with_replacement
from LammpsTreatmentFuncs import clean_data, clean_settings, find_sections, get_data, add_section_keyword, get_coeff, save_text_file

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

def file_unifier(directory, coeffsFile, dataList):
    # Class for handling Lammps data
    class Data:
        def __init__(self, data):
            self.data = data
            self.sectionIndexList = find_sections(self.data)

            # Header data
            self.header = self.data[:14]

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
            new_index_range = list(range(1, len(mass_types)+1))
            new_index_range = [str(val) for val in new_index_range]
            mass_zip = zip(mass_types, new_index_range)
            mass_change_dict = dict(mass_zip)
            
            # Change atom types to new types
            for massList in valid_masses:
                massList[0] = mass_change_dict[massList[0]]
            
            # Update self and add section keyword
            self.masses = add_section_keyword('Masses', valid_masses)

            return mass_change_dict
        
        def change_atom_types(self, mass_change_dict):
            for atomList in self.atoms:
                atomList[2] = mass_change_dict[atomList[2]]

            add_section_keyword('Atoms', self.atoms)

        def change_section_types(self, unioned_types, data_section):
            '''
            This function is different from change_mass_types as there is no removal of lines
            '''
            # Build new dict with old type keys and new type values
            new_index_range = list(range(1, len(unioned_types) + 1))
            new_index_range = [str(val) for val in new_index_range]
            type_zip = zip(unioned_types, new_index_range)
            type_change_dict = dict(type_zip)

            # Update data
            sectionData = getattr(self, data_section)
            for entryList in sectionData:
                entryList[1] = type_change_dict[entryList[1]]

            # Add section keywords
            sentenceCaseSection = data_section.capitalize()
            add_section_keyword(sentenceCaseSection, sectionData)

            return type_change_dict

        def change_header(self, typeList):
            self.header = [line.split() for line in self.header]
            for index, sectionType in enumerate(typeList):
                self.header[6+index][0] = sectionType

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

    def union_types(typeAttr, lammpsData=lammpsData):
        lammpsTypes = []
        for data in lammpsData:
            # Get attribute for determining types
            func = getattr(data, typeAttr)
            # Run function and append to list
            lammpsTypes.append(func())

        # Union sets to remove duplicates and sort into numerical order list 
        types = sorted(set().union(*lammpsTypes))
        numTypes = str(len(types))
        
        print(f'{typeAttr}\n Types: {types}\n Count: {numTypes}')

        return types, numTypes

    # Union sets and create sorted list for each type
    atomTypes, numAtomTypes = union_types('get_atom_types')
    bondTypes, numBondTypes = union_types('get_bond_types')
    angleTypes, numAngleTypes = union_types('get_angle_types')
    dihedralTypes, numDihedralTypes = union_types('get_dihedral_types')
    improperTypes, numImproperTypes = union_types('get_improper_types')

    # Update sections
    for data in lammpsData:
        bondDict = data.change_section_types(bondTypes, 'bonds')
        angleDict = data.change_section_types(angleTypes, 'angles')
        dihedralDict = data.change_section_types(dihedralTypes, 'dihedrals')
        improperDict = data.change_section_types(improperTypes, 'impropers')
        massDict = data.change_mass_types(atomTypes)
        data.change_atom_types(massDict)

    # Update header
    sectionTypeCounts = [numAtomTypes, numBondTypes, numAngleTypes, numDihedralTypes, numImproperTypes]
    for data in lammpsData:
        data.change_header(sectionTypeCounts)
        # Reinsert spaces between subsections of header
        data.header.insert(11, '\n')
        data.header.insert(6, '\n')
        data.header.insert(1, '\n')

    # Save data files
    
    for index, data in enumerate(lammpsData):
        # Combine all different data sources into one list
        combinedData = [data.header, data.masses, data.atoms, data.bonds, data.angles, data.dihedrals, data.impropers]
        # Flatten list of lists by one
        combinedData = [val for sublist in combinedData for val in sublist]

        # Save to text file
        save_text_file('cleaned' + dataList[index], combinedData)
    
    ####SETTINGS####

    # Load dataFile into python as a list of lists
    with open(coeffsFile, 'r') as f:
        settings = f.readlines()
    
    # Tidy settings and split
    settings = clean_settings(settings)
    settings = [line.split() for line in settings]

    # Create original atom type pair_coeff pairs
    originalPairTuples = list(combinations_with_replacement(atomTypes, 2))
    # Get all pair_coeffs
    pairCoeff = get_coeff("pair_coeff", settings)
    
    # Find valid pair_coeff pairs that are needed for this molecule
    validPairCoeff = []
    for pair in originalPairTuples:
        for coeff in pairCoeff:
            if coeff[1] == pair[0] and coeff[2] == pair[1]:
                validPairCoeff.append(coeff)
                break

    # Update atom types in pair_coeffs with massDict
    for pair in validPairCoeff:
        pair[1] = massDict[pair[1]]
        pair[2] = massDict[pair[2]]

    def valid_coeffs(coeffType, updateDict, settingsData=settings):
        # Get coeff lines
        coeffs = get_coeff(coeffType, settingsData)

        # Find valid coeffs from keys of updateDict
        validCoeffs = []
        for key in updateDict.keys():
            for coeff in coeffs:
                if coeff[1] == key:
                    validCoeffs.append(coeff)
                    break

        # Update coeffs with values of updateDict
        for coeff in validCoeffs:
            coeff[1] = updateDict[coeff[1]]

        return validCoeffs

    # Update coeff values
    validBondCoeff = valid_coeffs('bond_coeff', bondDict)
    validAngleCoeff = valid_coeffs('angle_coeff', angleDict)
    validDihedralCoeff = valid_coeffs('dihedral_coeff', dihedralDict)
    validImproperCoeff = valid_coeffs('improper_coeff', improperDict)

    # Combine all the coeff sources
    combinedCoeffs = [validPairCoeff, validBondCoeff, validAngleCoeff, validDihedralCoeff, validImproperCoeff]
    # Flatten list of lists by one
    combinedCoeffs = [val for sublist in combinedCoeffs for val in sublist]


    save_text_file('cleaned' + coeffsFile, combinedCoeffs)
