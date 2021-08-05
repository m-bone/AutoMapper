##############################################################################
# Developed by: Matthew Bone
# Last Updated: 30/07/2021
# Updated by: Matthew Bone
#
# Contact Details:
# Bristol Composites Institute (BCI)
# Department of Aerospace Engineering - University of Bristol
# Queen's Building - University Walk
# Bristol, BS8 1TR
# U.K.
# Email - matthew.bone@bristol.ac.uk
#
# File Description:
# Processes any number of LAMMPS 'read_data' input files, and a coefficients 
# file, in unison to reduce coefficient values down to the lowest possible 
# values. This was originally designed to work with Moltemplate system.data and
# system.in.settings files.
# If only one data file is specified this function is analogous to
# cleanup_moltemplate.sh.

# Assumptions:
# LAMMPS Atom Type is full
# Unit cell size is identical between files
# Pair coeffs are defined with numbers and not wildcarded with *
##############################################################################

import os
from natsort import natsorted
from itertools import combinations_with_replacement
from LammpsTreatmentFuncs import clean_data, clean_settings, add_section_keyword, save_text_file
from LammpsSearchFuncs import get_data, get_coeff, find_sections, get_header, convert_header

def file_unifier(directory, coeffsFile, dataList):
    # Go to file directory
    os.chdir(directory)

    # Load files from dataList, tidy and initialise class object
    lammpsData = []
    for dataFile in dataList:
        with open(dataFile, 'r') as f:
            data = f.readlines()

        # Tidy data
        data = clean_data(data)

        headerDict = get_header(data)

        # Initialise data class
        data = Data(data, headerDict)
        lammpsData.append(data)

    def union_types(typeAttr, lammpsData=lammpsData):
        lammpsTypes = []
        for data in lammpsData:
            # Get attribute for determining types
            getTypeAttr = 'get_' + typeAttr
            func = getattr(data, getTypeAttr)
            # Run function and append to list
            lammpsTypes.append(func())

        # Union sets to remove duplicates and sort into numerical order list 
        types = natsorted(set().union(*lammpsTypes))
        numTypes = (typeAttr, str(len(types))) # Tuple so that type can be accessed in dict later
        
        # Print number of types changed
        # print(f'{typeAttr}\n Types: {types}\n Count: {numTypes}') 
        
        return types, numTypes

    # Union sets and create sorted list for each type
    atomTypes, numAtomTypes = union_types('atom_types')
    bondTypes, numBondTypes = union_types('bond_types')
    angleTypes, numAngleTypes = union_types('angle_types')
    dihedralTypes, numDihedralTypes = union_types('dihedral_types')
    improperTypes, numImproperTypes = union_types('improper_types')

    # Update sections
    for data in lammpsData:
        bondDict = data.change_section_types(bondTypes, 'bonds')
        angleDict = data.change_section_types(angleTypes, 'angles')
        dihedralDict = data.change_section_types(dihedralTypes, 'dihedrals')
        improperDict = data.change_section_types(improperTypes, 'impropers')
        massDict = data.change_mass_types(atomTypes)
        data.change_atom_types(massDict)

    # Update header - will delete multiline comments and leave only the first
    sectionTypeCounts = [numAtomTypes, numBondTypes, numAngleTypes, numDihedralTypes, numImproperTypes]
    for data in lammpsData:
        data.change_header(sectionTypeCounts)

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
    # Currently, the h-bond flag value in hbond/dreiding is sorted as H_HB will always be 2 if H is in system
    # Apart from water or peroxide...
    validPairCoeff = []
    for pair in originalPairTuples:
        for coeff in pairCoeff:
            if coeff[1] == pair[0] and coeff[2] == pair[1]:
                validPairCoeff.append(coeff)

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

    # Save coeff file
    save_text_file('cleaned' + coeffsFile, combinedCoeffs)

# Class for handling Lammps data
class Data:
    def __init__(self, data, headerDict):
        self.data = data
        self.sectionIndexList = find_sections(self.data)

        # Header data
        self.header = headerDict

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
            # Add space to start of comment, if present, so it matches moltemplate
            if len(massList) >2:
                massList[2] = massList[2].rjust(2)
        
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
        # Iterate through type tuples and update header
        for typeData in typeList:
            self.header[typeData[0]] = [typeData[1]] # Must be list or >1 digit types get space separated. Str type shouldn't be a problem
        
        # Convert list values back to strings
        stringHeader = convert_header(self.header)
        
        # Restore header to a list of lists of strings
        self.header = stringHeader