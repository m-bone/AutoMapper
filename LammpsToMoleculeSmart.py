##############################################################################
# Developed by: Matthew Bone
# Last Updated: 08/02/2021
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
# Designed to work with LAMMPS data generated from Moltemplate, but can work
# for any LAMMPS data if it satisfies the header format in Assumptions.

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
import re
import sys
from LammpsTreatmentFuncs import clean_data, find_sections, get_data, add_section_keyword, save_text_file

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

    # Get first part of the LAMMPS header
    header = tidiedLines[1:6]
    header.insert(0, '\n')
    header =[[item] for item in header]

    # Get atoms data
    atoms = get_data('Atoms', tidiedLines, sectionIndexList)

    # Get bonds data
    bonds = get_data('Bonds', tidiedLines, sectionIndexList)
    
    # Search bond pair
    def pair_search(bond, bondAtom):
        if bond[2] == bondAtom:
            return bond[3]
        elif bond[3] == bondAtom:
            return bond[2]

    def search_loop(bonds, bondAtom, validAtomSet):
        nextBondAtomList = []

        for searchAtom in bondAtom:
            for bond in bonds:
                nextAtomID = pair_search(bond, searchAtom)
                if nextAtomID is not None:
                    validAtomSet.add(nextAtomID)
                    nextBondAtomList.append(nextAtomID)
                    # print(nextAtomID)
        
        return nextBondAtomList
            

    # Find bonds involving bonding atoms
    bondDistance = 3

    validAtomSet = set(bondingAtoms)
    for bondAtom in bondingAtoms:

        # Make bondAtom a list
        newBondAtomList = [bondAtom]
        
        i = 1
        while i <= bondDistance:
            newBondAtomList = search_loop(bonds, newBondAtomList, validAtomSet)
            i += 1

    # Need newBondAtomList to update with all the new atoms that need to be found and then search through bonds again for those numbers

    # validAtomSet = set(bondingAtoms)
    # for bondAtom in bondingAtoms:
    #     print(bondAtom)
    #     for bond in bonds:
    #         nextAtomID = pair_search(bond, bondAtom)
    #         if nextAtomID is not None:
    #             validAtomSet.add(nextAtomID)
    #             print(nextAtomID)
            

    checkList = sorted([int(val) for val in validAtomSet])
    print(checkList)
    
    
    # bonds = add_section_keyword('Bonds', bonds)

    # # Get angles data
    # angles = get_data('Angles', tidiedLines, sectionIndexList)
    # angles = add_section_keyword('Angles', angles)

    # # Get dihedrals
    # dihedrals = get_data('Dihedrals', tidiedLines, sectionIndexList)
    # dihedrals = add_section_keyword('Dihedrals', dihedrals)

    # # Get impropers
    # impropers = get_data('Impropers', tidiedLines, sectionIndexList)
    # impropers = add_section_keyword('Impropers', impropers)

    # # Rearrange atom data to get types, charges, coords - assume atom type full very important
    # types = [[atom[0], atom[2]] for atom in atoms]
    # types = add_section_keyword('Types', types)

    # charges = [[atom[0], atom[3]] for atom in atoms]
    # charges = add_section_keyword('Charges', charges)

    # coords = [[atom[0], atom[4], atom[5], atom[6]] for atom in atoms]
    # coords = add_section_keyword('Coords', coords)

    # # Combine to one long output list
    # outputList = []
    # totalList = [header, types, charges, coords, bonds, angles, dihedrals, impropers]
    
    # for keyword in totalList:
    #     outputList.extend(keyword)
        
    # # Output as text file
    # save_text_file(saveName + 'molecule.data', outputList)

lammps_to_molecule('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Pre_Reaction', ['1'], 'system.data', 'changed') # '33', '62'