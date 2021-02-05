import sys
from LammpsUnifiedCleaner import file_unifier
from LammpsToMolecule import lammps_to_molecule

directory = sys.argv[1]

# Unified data file clean
if sys.argv[2] == "clean":  
    coeffsFileName = sys.argv[3]
    dataFiles = sys.argv[4:]
    print(f'DataFiles List: {dataFiles}')

    file_unifier(directory, coeffsFileName, dataFiles)

# Produce molecule data file
elif sys.argv[2] == "molecule":
    cmdFileName = sys.argv[3]
    cmdSaveName = sys.argv[4]
    lammps_to_molecule(directory, cmdFileName , cmdSaveName)