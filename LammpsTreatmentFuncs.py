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
# A range of functions designed to pre-process, search, and save LAMMPS
# 'read_data' input files. Built for Bond_React processing but can be applied
# to a wide range of problems
##############################################################################

import re # For clean_data, clean_settings
from collections import Counter # For refine_data
from operator import itemgetter # For refine_data
from natsort import natsorted # For refine_data

# Function maybe moved to general function file later
def clean_data(lines):
    # Remove blank lines 
    lines = [line for line in lines if line != '\n']

    # Remove comments - negative lookbehind means label comments in masses are kept e.g # C_3
    lines = [re.sub(r'(?<!\d\s\s)#(.*)', '' ,line) for line in lines]

    # Remove 5 spaces moltemplate/LAMMPS puts in front of header file
    lines = [re.sub(r'\s{5}', '' ,line) for line in lines]

    # Remove newline terminators
    lines = [re.sub(r'\n', '', line) for line in lines]

    # Remove empty strings in list caused by comments being removed
    lines = [line for line in lines if line != '']

    # Remove trailing whitespace
    lines = [re.sub(r'\s+$', '', line) for line in lines]

    return lines

def clean_settings(lines):
    # Remove newline terminators
    lines = [re.sub(r'\n', '', line) for line in lines]

    # Remove tabs
    lines = [re.sub(r'\t', '', line) for line in lines]
    
    # Replace multiple whitespaces with one
    lines = [re.sub(r'\s{2,}', ' ', line) for line in lines]

    return lines

def find_sections(lines):
    # Find index of section keywords - isalpha works as no spaces, newlines or punc in section keywords
    sectionIndexList = [lines.index(line) for line in lines if line.isalpha()]

    # Add end of file as last index
    sectionIndexList.append(len(lines))

    return sectionIndexList

# Get data
def get_data(sectionName, lines, sectionIndexList):
    # Checks that section name is existing in LAMMPS data
    try:
        startIndex = lines.index(sectionName)
    except ValueError:
        # If doesn't exist, return empty list that can be added as normal to main list later
        data = []
        return data

    endIndex = sectionIndexList[sectionIndexList.index(startIndex) + 1]
    
    data = lines[startIndex+1:endIndex] # +1 means sectionName doesn't get included
    data = [val.split() for val in data]

    return data

def refine_data(data, searchIndex: list, IDset):
    '''
    Search multiple indices for matching atomID value.
    If match is found keep that list row in the data.
    Only output that row if the row appears len(searchIndex) times
    This means the data row contains valid atomIDs in all possible ID positions
    '''
    # Convenience to convert int values to list
    if type(searchIndex) is not list:
        searchIndex = [searchIndex]

    # Adds rows that include a valid atom ID to list
    possibleValidData = []
    for atomID in IDset:
        for row in data:
            for index in searchIndex:
                if row[index] == atomID:
                    possibleValidData.append(row)

    # Lammps IDs found in above search
    possibleValidIDs = [val[0] for val in possibleValidData]
    IDCount = dict(Counter(possibleValidIDs))
    # If ID counter is the same size as the search index, ID is valid and gets added to data
    validIDs = [key for key in IDCount.keys() if IDCount[key] == len(searchIndex)]        
    validData = []
    for row in possibleValidData:
        if row[0] in validIDs:
            validData.append(row)
            # Stops duplicate IDs being refound in the future
            validIDIndex = validIDs.index(row[0])
            del validIDs[validIDIndex]

    # Re-sort validData by ID, use natsort as values are str not int
    validData = natsorted(validData, key=itemgetter(0))
    
    return validData

def get_coeff(coeffName, settingsData):
    # Inputs pre-split data
    # Return all lines that include coeffName in the [0] index
    coeffs = [line for line in settingsData if line[0] == coeffName]
    
    return coeffs

def add_section_keyword(sectionName, data):
    # Don't add keyword if list is empty - empty list means no section in file
    if len(data) == 0:
        return data

    # Add keyword name to start of list
    data.insert(0, '\n')
    data.insert(1, [sectionName])
    data.insert(2, '\n')

    return data

def save_text_file(fileName, dataSource):
    # Save to text file
    with open(fileName, 'w') as f:
        for item in dataSource:
            line = " ".join(item)
            if line != '\n':
                f.write("%s\n" % line)
            else:
                f.write(line)

# Create comment string with bond atoms and edge atoms
def format_comment(IDlist, comment):
    atomList = natsorted(list(IDlist))
    atomList.insert(0, comment)
    atomString = [' '.join(atomList)] # Has to be list of lists to pass through later code

    return atomString

# Search bond pair
def pair_search(bond, bondAtom):
    if bond[2] == bondAtom:
        return bond[3]
    elif bond[3] == bondAtom:
        return bond[2]

# Loop through atomIDs, possible bonds and find valid bonds
def search_loop(bonds, bondAtom):
    nextBondAtomList = []

    for searchAtom in bondAtom:
        for bond in bonds:
            nextAtomID = pair_search(bond, searchAtom)
            if nextAtomID is not None:
                nextBondAtomList.append(nextAtomID)
    
    return nextBondAtomList
        
def find_partial_structure(bondingAtoms, originalBonds, bondDistance=3):
# Find bonds involving bonding atoms
    validAtomSet = set(bondingAtoms)
    edgeAtomSet = set()
    for bondAtom in bondingAtoms:

        # Make bondAtom a list
        newBondAtomList = [bondAtom]
        
        i = 1
        while i <= bondDistance:
            newBondAtomList = search_loop(originalBonds, newBondAtomList)
            if i == 1:
                # Stop search from finding other bonding atom if they are bound together
                newBondAtomList = [val for val in newBondAtomList if val not in bondingAtoms]

            if i < bondDistance:
                # Add list as individual elements
                [validAtomSet.add(val) for val in newBondAtomList]

            else:
                # Determine which of the last obtained atom IDs have more bonds
                # These atoms should be edge atoms
                furthestAtomSet = set(newBondAtomList)
                possibleEdgeAtoms = furthestAtomSet.difference(validAtomSet)
                
                [validAtomSet.add(val) for val in newBondAtomList]

                # Run another loop to determine if possibleEdgeAtoms have other bonds
                for searchAtom in possibleEdgeAtoms:
                    bondCount = 0
                    for bond in originalBonds:
                        nextAtomID = pair_search(bond, searchAtom)
                        if nextAtomID is not None:
                            bondCount += 1
                    if bondCount > 1: # All atoms will have at least one bond
                        edgeAtomSet.add(searchAtom)
            
            # Increment iterator
            i += 1
    return validAtomSet, edgeAtomSet