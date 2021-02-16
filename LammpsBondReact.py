##############################################################################
# Developed by: Matthew Bone
# Last Updated: 16/02/2021
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
# This script is designed to be run from run_bond_react.sh, allowing the user
# to call Bond_React python code from wherever their data is stored
##############################################################################

import sys
from LammpsUnifiedCleaner import file_unifier
from LammpsToMolecule import lammps_to_molecule
from LammpsToMoleculePartial import lammps_to_molecule_partial
from LammpsToLammpsPartial import lammps_to_lammps_partial

directory = sys.argv[1]

# Unified data file clean
if sys.argv[2] == "clean":  
    coeffsFileName = sys.argv[3]
    dataFiles = sys.argv[4:]
    print(f'DataFiles List: {dataFiles}')

    file_unifier(directory, coeffsFileName, dataFiles)

# Produce molecule data file
elif sys.argv[2] == "molecule":
    dataFileName = sys.argv[3]
    molSaveName = sys.argv[4]
    lammps_to_molecule(directory, dataFileName , molSaveName)

# Produce partial molecule data file
elif sys.argv[2] == "molecule-partial":
    dataFileName = sys.argv[3]
    molSaveName = sys.argv[4]
    bondingAtoms = sys.argv[5:]
    lammps_to_molecule_partial(directory, dataFileName, molSaveName, bondingAtoms)

# Produce partial lammps data file
elif sys.argv[2] == "lammps-partial":
    dataFileName = sys.argv[3]
    molSaveName = sys.argv[4]
    bondingAtoms = sys.argv[5:]
    lammps_to_lammps_partial(directory, dataFileName, molSaveName, bondingAtoms)