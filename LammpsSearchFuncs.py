##############################################################################
# Developed by: Matthew Bone
# Last Updated: 23/02/2021
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
# A range of functions designed to search LAMMPS files for information.
# These functions work for a 'read_data' files and 'molecule' files
##############################################################################
import re

# Get data
def get_data(sectionName, lines, sectionIndexList):
    # Checks that section name is existing in LAMMPS data
    startIndex = lines.index(sectionName)

    endIndex = sectionIndexList[sectionIndexList.index(startIndex) + 1]
    
    data = lines[startIndex+1:endIndex] # +1 means sectionName doesn't get included
    data = [val.split() for val in data]

    return data

def get_coeff(coeffName, settingsData):
    # Inputs pre-split data
    # Return all lines that include coeffName in the [0] index
    coeffs = [line for line in settingsData if line[0] == coeffName]
    
    return coeffs

def get_top_comments(data):
    '''
    Finds the top comments in partial data files and does some cleaning
    These comments are used in other automation processes
    '''

    topComments = []

    for line in data:
        if '#' in line:
            # Remove \n, the # and split line by spaces
            cleanedLine = re.sub(r'\n', '', line)
            cleanedLine = re.sub(r'#', '', cleanedLine)
            cleanedLine = cleanedLine.split()

            topComments.append(cleanedLine)
        else:
            # Stops it finding comments later in the file
            break

    return topComments

def find_sections(lines):
    # Find index of section keywords - isalpha works as no spaces, newlines or punc in section keywords
    sectionIndexList = [lines.index(line) for line in lines if line.isalpha()]

    # Add end of file as last index
    sectionIndexList.append(len(lines))

    return sectionIndexList

# Search bond pair
def pair_search(bond, bondAtom):
    '''
    Check if either atomID in a bond is the desired atomID.
    Will return None if no match is found.
    '''
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
    edgeAtomList = []

    for bondAtom in bondingAtoms:

        # Make bondAtom a list
        newBondAtomList = [bondAtom]
        
        i = 1
        while i <= bondDistance:
            newBondAtomList = search_loop(originalBonds, newBondAtomList)
            if i == 1: # First pass - Stop search from finding other bonding atom if they are bound together
                newBondAtomList = [val for val in newBondAtomList if val not in bondingAtoms]

            if i < bondDistance: # Before bond distance is reached
                # Add list as individual elements 
                for atom in newBondAtomList:
                    validAtomSet.add(atom)

            else: # Once bond distance is reached
                # Determine which of the last obtained atom IDs have further bonds               
                # newBondAtomList at this point contains edge atoms of an order, and other atoms found before
                possibleEdgeAtoms = [val for val in newBondAtomList if val not in validAtomSet]

                # Add list as individual elements - has to be after possibleEdgeAtoms
                for atom in newBondAtomList:
                    validAtomSet.add(atom)

                # Run another loop to determine if possibleEdgeAtoms have other bonds
                for searchAtom in possibleEdgeAtoms:
                    bondCount = 0
                    for bond in originalBonds:
                        nextAtomID = pair_search(bond, searchAtom)
                        if nextAtomID is not None:
                            bondCount += 1
                    if bondCount > 1: # All atoms will have at least one bond
                        edgeAtomList.append(searchAtom)
            
            # Increment iterator
            i += 1

    # Get edge atom neighbours
    edgeAtomFingerprintDict = get_neighbours(edgeAtomList, originalBonds)
    
    # Filter out validAtomIDs that are within the partial structure
    filteredFingerprintDict = {}
    for key, atomList in edgeAtomFingerprintDict.items():
        cutList = [atom for atom in atomList if atom not in validAtomSet]
        filteredFingerprintDict[key] = cutList

    return validAtomSet, edgeAtomList, filteredFingerprintDict

def get_neighbours(atomIDList, bondsList):
    boundAtomsList = []

    # Determine what atoms are bound to an initial atom
    for atom in atomIDList:
        bondingAtoms = []
        for bond in bondsList:
            pairResult = pair_search(bond, atom)
            if pairResult is not None:
                bondingAtoms.append(pairResult)

        boundAtomsList.append([atom, bondingAtoms])

    # Create dictionary of initial atom keys and bound atom list values
    boundAtomsDict = {val[0]: val[1] for val in boundAtomsList}

    return boundAtomsDict