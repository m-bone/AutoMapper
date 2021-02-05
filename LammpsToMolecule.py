import os
import re
import sys
from LammpsTreatmentFuncs import clean_data, find_sections, get_data, add_section_keyword, save_text_file


# ASSUMPTIONS
# LAMMPS Atom Type is full

def lammps_to_molecule(directory, fileName, saveName):
    # Go to file directory
    os.chdir(directory)

    # Load file into python as a list of lists
    with open(fileName, 'r') as f:
        lines = f.readlines()
    

    # Tidy input
    tidiedLines = clean_data(lines)
    
    # Build sectionIndexList
    sectionIndexList = find_sections(tidiedLines)

    # Get first part of the LAMMPS header
    header = tidiedLines[1:6]
    header.insert(0, '\n')
    header =[[item] for item in header]

    # Get atoms data
    atoms = get_data('Atoms', tidiedLines, sectionIndexList)

    # Get bonds data
    bonds = get_data('Bonds', tidiedLines, sectionIndexList)
    bonds = add_section_keyword('Bonds', bonds)

    # Get angles data
    angles = get_data('Angles', tidiedLines, sectionIndexList)
    angles = add_section_keyword('Angles', angles)

    # Get dihedrals
    dihedrals = get_data('Dihedrals', tidiedLines, sectionIndexList)
    dihedrals = add_section_keyword('Dihedrals', dihedrals)

    # Get impropers
    impropers = get_data('Impropers', tidiedLines, sectionIndexList)
    impropers = add_section_keyword('Impropers', impropers)

    # Rearrange atom data to get types, charges, coords - assume atom type full very important
    types = [[atom[0], atom[2]] for atom in atoms]
    types = add_section_keyword('Types', types)

    charges = [[atom[0], atom[3]] for atom in atoms]
    charges = add_section_keyword('Charges', charges)

    coords = [[atom[0], atom[4], atom[5], atom[6]] for atom in atoms]
    coords = add_section_keyword('Coords', coords)

    # Combine to one long output list
    outputList = []
    totalList = [header, types, charges, coords, bonds, angles, dihedrals, impropers]
    
    for keyword in totalList:
        outputList.extend(keyword)
        
    # Output as text file
    save_text_file(saveName + 'molecule.data', outputList)

