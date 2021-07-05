#!/usr/bin/env python3
##############################################################################
# Developed by: Matthew Bone
# Last Updated: 16/06/2021
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

import argparse
from LammpsUnifiedCleaner import file_unifier
from LammpsToMolecule import lammps_to_molecule
from LammpsToMoleculePartial import lammps_to_molecule_partial
from LammpsToLammpsPartial import lammps_to_lammps_partial
from PathSearch import map_from_path
from MapProcessor import map_processor, save_output

# Init the parser
parser = argparse.ArgumentParser(description='Run preprocessing tools for LAMMPS simulation using fix bond/react')

# List of arguments for command line
parser.add_argument('directory', metavar='directory', type=str, nargs=1, help='Directory of file(s), can be found in bash with . or $PWD')
parser.add_argument('tool', metavar='tool', type=str, nargs=1, choices=['clean', 'molecule', 'molecule-partial', 'lammps-partial', 'map', 'map-only', 'test'], help='Name of tool to be used. Possible tools: clean, molecule, molecule-partial, lammps-partial, map, map-only')
parser.add_argument('data_files', metavar='data_files', nargs='+', help='Name of file(s) to be acted on. If tool is "map" then this must be two files ordered as pre-bond post-bond. If "clean" this can be a list of files in any order')
parser.add_argument('--coeff_file', metavar='coeff_file', nargs=1, help='Argument for the "clean" tool: a coefficients file to be cleaned')
parser.add_argument('--save_name', metavar='save_name', nargs='+', help='Argument for "molecule", "molecule-partial" and "lammps-partial" tools: a prefix that is added to "molecule.data" for the file name of the new file')
parser.add_argument('--ba', metavar='bonding_atoms', nargs='+', help='Argument for "molecule", "molecule-partial", "lammps-partial" and "map" tools: atom IDs of the atoms that will be involved in creating a new bond, separated by white space. Order of atoms must be the same between molecule files, when mapping.')
parser.add_argument('--ebt', metavar='elements_by_type', nargs='+', help='Argument for "molecule-partial", "lammps-partial", "map" and "map-only" tools: series of elements symbols in the same order as the types specified in the data file and separated with a whitespace')
parser.add_argument('--da', metavar='delete_atoms', nargs='+', help='An optional argument for "molecule", "molecule-partial", "lammps-partial" and "map" tools: atom IDs of the atoms that will be deleted after the bond has formed, separated by white space')
parser.add_argument('--debug', action='store_true', help='Argument for "map" and "map-only": prints debugging statements highlighting which parts of the path search determine a mapped atom pair')

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

if tool != 'clean' and tool != 'map' and len(args.data_files) > 1:
    parser.error('This tool can only take 1 data_file as input')

if 'map' in tool and (len(args.data_files) != 2 or args.ebt is None):
    parser.error('Map tools require 2 data_files and --ebt (elements by type) arguments')

if tool == 'map' and (len(args.save_name) != 2 or len(args.ba) != 4):
    parser.error('Map tool requires 2 save_names (for pre and post molecule files) and --ba (bonding atoms) with 4 atomIDs specified')

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

# Combined molecule-partial and map creation code - can handle edge atom moving
elif tool == 'map':
    map_processor(directory, args.data_files[0], args.data_files[1], args.save_name[0], args.save_name[1], args.ba[:2], args.ba[2:], args.da, args.ebt, args.debug)

# Produce just a map file
elif tool == 'map-only':
    mappedIDList, _, commentPreBondingAtoms, commentPreEdgeAtomDict, commentPreDeleteAtoms = map_from_path(directory, args.data_files[0], args.data_files[1], args.ebt, args.debug)
    save_output(mappedIDList, commentPreBondingAtoms, commentPreEdgeAtomDict, commentPreDeleteAtoms)
