#!/usr/bin/env python3
##############################################################################
# Developed by: Matthew Bone
# Last Updated: 09/04/2021
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

import argparse
from LammpsUnifiedCleaner import file_unifier
from LammpsToMolecule import lammps_to_molecule
from LammpsToMoleculePartial import lammps_to_molecule_partial
from LammpsToLammpsPartial import lammps_to_lammps_partial

# Init the parser
parser = argparse.ArgumentParser(description='Run preprocessing tools for LAMMPS simulation using fix bond/react')

# List of arguments for command line
parser.add_argument('directory', metavar='directory', type=str, nargs=1, help='Directory of file(s), can be found in bash with . or $PWD')
parser.add_argument('tool', metavar='tool', type=str, nargs=1, choices=['clean', 'molecule', 'molecule-partial', 'lammps-partial'], help='Name of tool to be used. Possible tools: clean, molecule, molecule-partial, lammps-partial')
parser.add_argument('data_files', metavar='data_files', nargs='+', help='Name of file to be acted on. If tool is "clean" this can be a list of files')
parser.add_argument('--coeff_file', metavar='coeff_file', nargs=1, help='Argument for the "clean" tool: a coefficients file to be cleaned')
parser.add_argument('--save_name', metavar='save_name', nargs=1, help='Argument for "molecule", "molecule-partial" and "lammps-partial" tools: sets the file name for the new file')
parser.add_argument('--ba', metavar='bonding_atoms', nargs=2, help='Argument for "molecule", "molecule-partial" and "lammps-partial" tools: atom IDs of the atoms that will be involved in creating a new bond')
parser.add_argument('--ebt', metavar='elements_by_type', nargs='+', help='Argument for "molecule-partial" and "lammps-partial" tools: list of elements symbols in the same order as the types specified in the data file')
parser.add_argument('--da', metavar='delete_atoms', nargs='+', help='Argument for "molecule", "molecule-partial" and "lammps-partial" tools: atom IDs of the atoms that will be deleted after the bond has formed')

# Get arguments from parser
args = parser.parse_args()
# Take compulsory args out of list - other args done later
tool = args.tool[0]
directory = args.directory[0]


# Throw errors if arguments are missing from certain tools
if tool == 'clean' and (args.coeff_file is None):
    parser.error('"clean" tool requries --coeff_file')

if tool == 'molecule' and (args.save_name is None or args.ba is None):
    parser.error('"molecule" tool requries --save_name and --ba (bonding atoms) arguments')

if 'partial' in tool and (args.save_name is None or args.ba is None or args.ebt is None):
    parser.error('"molecule-partial" or "lammps-partial" tools require --save_name, --ba (bonding atoms) and --ebt (elements by type) arguments')

if tool != 'clean' and len(args.data_files) > 1:
    parser.error('Only the "clean" tool can take more than 1 data_file as input')


# Unified data file clean
if tool == "clean":  
    print(f'DataFiles List: {args.data_files}')
    file_unifier(directory, args.coeff_file[0], args.data_files)

# Produce molecule data file
elif tool == "molecule":
    lammps_to_molecule(directory, args.data_files[0], args.save_name[0], args.ba, args.da)

# Produce partial molecule data file
elif tool == "molecule-partial":
    lammps_to_molecule_partial(directory, args.data_files[0], args.save_name[0], args.ebt, args.ba, args.da)

# Produce partial lammps data file
elif tool == "lammps-partial":
    lammps_to_lammps_partial(directory, args.data_files[0], args.save_name[0], args.ebt, args.ba, args.da)