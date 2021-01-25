import os
import re
import sys
from LammpsTreatmentFuncs import clean_data, clean_settings, find_sections, get_data, add_section_keyword

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

def first_file_unifier(directory, dataFile1, dataFile2, settingsFile):
    # Go to file directory
    os.chdir(directory)

    # Load dataFile into python as a list of lists
    with open(dataFile1, 'r') as f:
        data1 = f.readlines()

    # Tidy data1
    data1 = clean_data(data1)


    # Load dataFile into python as a list of lists
    with open(dataFile2, 'r') as f:
        data2 = f.readlines()

    # Tidy data2
    data2 = clean_data(data2)


    # Load dataFile into python as a list of lists
    with open(settingsFile, 'r') as f:
        settings = f.readlines()
    
    # Tidy settings
    settings = clean_settings(settings)

    # Class for handling Lammps data
    class Data:
        def __init__(self, data):
            self.data = data
            self.sectionIndexList = find_sections(self.data)
            self.atoms = get_data('Atoms', self.data, self.sectionIndexList)
            self.masses = get_data('Masses', self.data, self.sectionIndexList)
        
        def get_atom_types(self):
            atom_types = {atom[2] for atom in self.atoms}
            return atom_types
        
        def get_mass_types(self, unioned_atom_types):
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
        
        def update_atoms(self, mass_change_dict):
            for atomList in self.atoms:
                atomList[2] = mass_change_dict[atomList[2]]

    # Initialise classes for Lammps data
    Lammps1 = Data(data1)
    Lammps2 = Data(data2)

    # Union sets and create sorted list
    atomTypes = Lammps1.get_atom_types().union(Lammps2.get_atom_types())
    atomTypes = sorted(atomTypes)

    # Build new masses section - eliminate unused masses
    processedMass, massDict = Lammps1.get_mass_types(atomTypes)
    processedMass = add_section_keyword('Masses', processedMass)

    # Update atom types and get atoms sections
    Lammps1.update_atoms(massDict)
    Lammps2.update_atoms(massDict)
    processedAtoms1 = add_section_keyword('Atoms', Lammps1.atoms)
    processedAtoms2 = add_section_keyword('Atoms', Lammps2.atoms)

    print()

first_file_unifier('/home/matt/Documents/Oct20-Dec20/Bonding_Test/React_Test/PythonCode', 'pre-system.data', 'post-system.data', 'system.in.settings')