##############################################################################
# Developed by: Matthew Bone
# Last Updated: 12/02/2021
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
# Converts LAMMPS 'read_data' input files into a cut down partial structure,
# still in LAMMPS data input file format. This is more for testing as these
# files can be read by visualisors like Ovito 

# Assumptions:
# LAMMPS Atom Type is full
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
##############################################################################

import os
from natsort import natsorted
from LammpsTreatmentFuncs import clean_data, find_sections, get_data, add_section_keyword, save_text_file, refine_data, find_partial_structure

def lammps_to_molecule(directory, bondingAtoms, fileName, saveName):
    # Go to file directory
    os.chdir(directory)

    # Load file into python as a list of lists
    with open(fileName, 'r') as f:
        lines = f.readlines()
    
    # Tidy input
    tidiedLines = clean_data(lines)
    
    # Build sectionIndexList
    sectionIndexList = find_sections(tidiedLines)

    # Get original bonds data
    originalBonds = get_data('Bonds', tidiedLines, sectionIndexList)
    
    validAtomSet, edgeAtomSet = find_partial_structure(bondingAtoms, originalBonds)
    
    # Get masses data
    masses = get_data('Masses', tidiedLines, sectionIndexList)
    masses = add_section_keyword('Masses', masses)

    # Get atoms data
    atoms = get_data('Atoms', tidiedLines, sectionIndexList)
    atoms = refine_data(atoms, 0, validAtomSet)
    atoms = add_section_keyword('Atoms', atoms)
    # SAYS IT HAS 12 ATOMS, ONLY OUTPUTS 9 - GETS COUNT WRONG ALL OVER THE PLACE
    
    # Get new bonds data
    bonds = refine_data(originalBonds, [2, 3], validAtomSet)
    bonds = add_section_keyword('Bonds', bonds)
    
    # Get angles data
    angles = get_data('Angles', tidiedLines, sectionIndexList)
    angles = refine_data(angles, [2, 3, 4], validAtomSet)
    angles = add_section_keyword('Angles', angles)

    # Get dihedrals
    dihedrals = get_data('Dihedrals', tidiedLines, sectionIndexList)
    dihedrals = refine_data(dihedrals, [2, 3, 4, 5], validAtomSet)
    dihedrals = add_section_keyword('Dihedrals', dihedrals)

    # Get impropers
    impropers = get_data('Impropers', tidiedLines, sectionIndexList)
    impropers = refine_data(impropers, [2, 3, 4, 5], validAtomSet)
    impropers = add_section_keyword('Impropers', impropers)

    # Get and change header values
    header = tidiedLines[1:14]
    header = [val.split() for val in header]
    # Update numbers with new lengths of data
    for index, data in enumerate([atoms, bonds, angles, dihedrals, impropers]):
        if len(data) > 0: # This corrects for the added section keyword lines
            header[index][0] = str(len(data) - 3)
        else:
            header[index][0] = str(0)
    header.insert(10, '\n')
    header.insert(5, '\n')
    header.insert(0, '\n')

    # Add edge atoms as comment in header
    edgeAtomList = natsorted(list(edgeAtomSet))
    edgeAtomList.insert(0, '#Edge Atoms:')
    edgeAtomString = [[' '.join(edgeAtomList)]] # Has to be list of lists to pass through later code
    
    # Combine to one long output list
    outputList = []
    totalList = [edgeAtomString, header, masses, atoms, bonds, angles, dihedrals, impropers]
    
    for keyword in totalList:
        outputList.extend(keyword)
        
    # Output as text file
    save_text_file(saveName + 'molecule.data', outputList)

lammps_to_molecule('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Pre_Reaction', ['28'], 'system.data', 'lammpstest')