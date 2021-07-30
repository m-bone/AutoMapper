#!/usr/bin/env python3
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
# This wrapper is used to call all the AutoMapper family of tools from a command
# line. Linux users can consider adding the location of this repo to their PATH
# so that this wrapper can be called from anywhere.
##############################################################################

import sys
# Check python version
pythonVersion = sys.version_info
if pythonVersion[0] == 2: # Python 2 check - might be needless due to the shebang
    print('This script is not compatible with Python 2. Please update to Python 3')
    sys.exit()

if pythonVersion[0] == 3 and pythonVersion[1] < 8: # Note on earlier versions of Python 3
    print('Note: This code was built in Python 3.8. It should work with earlier Python versions, but performance may vary.')

# Check that natsort is installed
try:
    import natsort
except ModuleNotFoundError:
    print('The package natsort was not found in the Python modules. Please install natsort before continuing.')
    sys.exit()

import argparse
from LammpsUnifiedCleaner import file_unifier
from LammpsToMolecule import lammps_to_molecule
from MapProcessor import map_processor

# Init the parser
parser = argparse.ArgumentParser(description='Run preprocessing tools for LAMMPS simulation using fix bond/react')

# List of arguments for command line
parser.add_argument('directory', metavar='directory', type=str, nargs=1, help='Directory of file(s), can be found in bash with . or $PWD')
parser.add_argument('tool', metavar='tool', type=str, nargs=1, choices=['clean', 'molecule', 'map'], help='Name of tool to be used. Possible tools: clean, molecule, map')
parser.add_argument('data_files', metavar='data_files', nargs='+', help='Name of file(s) to be acted on. If tool is "map" then this must be two files ordered as pre-bond post-bond. If "clean" this can be a list of files in any order')
parser.add_argument('--coeff_file', metavar='coeff_file', nargs=1, help='Argument for the "clean" tool: a coefficients file to be cleaned')
parser.add_argument('--save_name', metavar='save_name', nargs='+', help='Argument for "molecule" and "map" tools: the file name of the new file(s)')
parser.add_argument('--ba', metavar='bonding_atoms', nargs='+', help='Argument for "molecule" and "map" tools: atom IDs of the atoms that will be involved in creating a new bond, separated by a space. Order of atoms must be the same between molecule files when mapping.')
parser.add_argument('--ebt', metavar='elements_by_type', nargs='+', help='Argument for the "map" tools: series of elements symbols in the same order as the types specified in the data file and separated with a space')
parser.add_argument('--da', metavar='delete_atoms', nargs='+', help='An optional argument for "molecule" and "map" tools: atom IDs of the atoms that will be deleted after the bond has formed, separated by a space')
parser.add_argument('--debug', action='store_true', help='An optional argument for the "map" tool: prints debugging statements with information on the path search and map processor.')

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

if tool == 'molecule' and len(args.data_files) > 1:
    parser.error('The molecule tool can only take 1 data_file as input')

if tool == 'map' and (len(args.data_files) != 2 or len(args.save_name) != 2):
    parser.error('The map tool requires 2 data_files and 2 save_names (for pre and post molecule files)')

if tool == 'map' and (len(args.ba) != 4 or args.ebt is None):
    parser.error('The map tool requires --ba (bonding atoms) with 4 atomIDs specified and --ebt (elements by type) arguments')

# Unified data file clean
if tool == "clean":  
    print(f'DataFiles List: {args.data_files}')
    file_unifier(directory, args.coeff_file[0], args.data_files)

# Produce molecule data file
elif tool == "molecule":
    lammps_to_molecule(directory, args.data_files[0], args.save_name[0], args.ba, args.da)

# Combined molecule and map creation code
elif tool == 'map':
    map_processor(directory, args.data_files[0], args.data_files[1], args.save_name[0], args.save_name[1], args.ba[:2], args.ba[2:], args.da, args.ebt, args.debug)

