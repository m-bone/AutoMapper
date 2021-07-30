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
# Converts LAMMPS 'read_data' input files into LAMMPS molecule format files.
# Should read any valid format of LAMMPS input file; header assumptions have
# been removed.
##############################################################################

import os
from LammpsTreatmentFuncs import clean_data, add_section_keyword, refine_data, save_text_file, format_comment
from LammpsSearchFuncs import get_data, find_sections, get_header, convert_header

def lammps_to_molecule(directory, fileName, saveName, bondingAtoms: list, deleteAtoms=None, validIDSet=None, renumberedAtomDict=None):
    # Go to file directory
    os.chdir(directory)

    # Load file into python as a list of lists
    with open(fileName, 'r') as f:
        lines = f.readlines()
    
    # Tidy input
    tidiedLines = clean_data(lines)
    
    # Build sectionIndexList
    sectionIndexList = find_sections(tidiedLines)

    # Get atoms data
    atoms = get_data('Atoms', tidiedLines, sectionIndexList)
    atoms = refine_data(atoms, 0, validIDSet, renumberedAtomDict)

    # Get bonds data
    bonds = get_data('Bonds', tidiedLines, sectionIndexList)
    bonds = refine_data(bonds, [2, 3], validIDSet, renumberedAtomDict)
    bondInfo = ('bonds', len(bonds))
    bonds = add_section_keyword('Bonds', bonds)

    # Get angles data
    angles = get_data('Angles', tidiedLines, sectionIndexList)
    angles = refine_data(angles, [2, 3, 4], validIDSet, renumberedAtomDict)
    angleInfo = ('angles', len(angles))
    angles = add_section_keyword('Angles', angles)

    # Get dihedrals
    dihedrals = get_data('Dihedrals', tidiedLines, sectionIndexList)
    dihedrals = refine_data(dihedrals, [2, 3, 4, 5], validIDSet, renumberedAtomDict)
    dihedralInfo = ('dihedrals', len(dihedrals))
    dihedrals = add_section_keyword('Dihedrals', dihedrals)

    # Get impropers
    impropers = get_data('Impropers', tidiedLines, sectionIndexList)
    impropers = refine_data(impropers, [2, 3, 4, 5], validIDSet, renumberedAtomDict)
    improperInfo = ('impropers', len(impropers))
    impropers = add_section_keyword('Impropers', impropers)

    # Rearrange atom data to get types, charges, coords - assume atom type full very important
    types = [[atom[0], atom[2]] for atom in atoms]
    typeInfo = ('atoms', len(types))
    types = add_section_keyword('Types', types)

    charges = [[atom[0], atom[3]] for atom in atoms]
    charges = add_section_keyword('Charges', charges)

    coords = [[atom[0], atom[4], atom[5], atom[6]] for atom in atoms]
    coords = add_section_keyword('Coords', coords)

    # Get and change header values
    header = get_header(tidiedLines)
    
    # Update numbers with new lengths of data if new IDs have been supplied
    if validIDSet is not None:
        for info in [typeInfo, bondInfo, angleInfo, dihedralInfo, improperInfo]:
            header[info[0]] = [info[1]]

    # Create bonding atom comment
    commentString = format_comment(bondingAtoms, 'Bonding_Atoms ')
    if deleteAtoms is not None:
        deleteAtomComment = format_comment(deleteAtoms, 'Delete_Atoms')
        commentString.extend([deleteAtomComment])
    header['comment'].extend(commentString)

    # Remove unnecessary header elements
    keepList = ['comment', 'atoms', 'bonds', 'angles', 'dihedrals', 'impropers']
    cutHeader = {key: header[key] for key in keepList}

    # Convert header back to list of lists of strings
    header = convert_header(cutHeader)

    # Combine to one long output list
    outputList = []
    totalList = [header, types, charges, coords, bonds, angles, dihedrals, impropers]
    
    for keyword in totalList:
        outputList.extend(keyword)
        
    # Output as text file
    save_text_file(saveName, outputList)

