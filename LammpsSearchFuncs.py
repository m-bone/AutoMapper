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
# A range of functions designed to search LAMMPS files for information.
# These functions work for 'read_data' files and 'molecule' files
##############################################################################
import re
from collections import Counter
from LammpsTreatmentFuncs import clean_data

# Get data
def get_data(sectionName, lines, sectionIndexList, useExcept = True):
    if useExcept: # Checks that section name is existing in LAMMPS data
        try:
            startIndex = lines.index(sectionName)
        except ValueError:
            # If doesn't exist, return empty list that can be added as normal to main list later
            data = []
            return data

    else: # Allows for later try/except blocks to catch missing section names
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

def read_top_comments(topComments, commentTerm):
    # Iterate through comments to find the one with the matching term
    # Return remaining list of values
    # If value is not in comments, None will be returned
    for comment in topComments:
        if commentTerm in comment:
            # Delete comment term from list
            cutComment = comment.copy()
            cutComment.remove(commentTerm)
            return cutComment

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
        
def find_partial_structure(bondingAtoms, originalBondList, deleteAtoms, bondDistance=3):
    # Find bonds within a specified distance of the bonding atoms
    
    # Convert bondingAtoms to list if string given
    if type(bondingAtoms) == str:
        bondingAtoms = [bondingAtoms]

    # Add delete atoms to valid atoms if present
    initialValidAtoms = bondingAtoms.copy()
    if deleteAtoms is not None:
        initialValidAtoms.extend(deleteAtoms) # Allows partial structure tools to work when byproducts are formed and deleted

    validAtomSet = set(initialValidAtoms)
    edgeAtomList = []

    for bondAtom in bondingAtoms:

        # Make bondAtom a list
        newBondAtomList = [bondAtom]
        
        i = 1
        while i <= bondDistance:
            newBondAtomList = search_loop(originalBondList, newBondAtomList)
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
                    for bond in originalBondList:
                        nextAtomID = pair_search(bond, searchAtom)
                        if nextAtomID is not None:
                            bondCount += 1
                    if bondCount > 1: # All atoms will have at least one bond
                        edgeAtomList.append(searchAtom)
            
            # Increment iterator
            i += 1

    return validAtomSet, edgeAtomList

def edge_atom_fingerprint_ids(edgeAtomList, originalBondList, validAtomSet):
    # Get edge atom neighbours
    edgeAtomFingerprintDict = get_neighbours(edgeAtomList, originalBondList, []) # Bonding atoms given as blank list, edge atoms can never have bonding atoms as a neighbour so not a problem
    
    # Filter out validAtomIDs that are within the partial structure
    filteredFingerprintDict = {}
    for key, atomList in edgeAtomFingerprintDict.items():
        cutList = [atom for atom in atomList if atom not in validAtomSet]
        filteredFingerprintDict[key] = cutList

    return filteredFingerprintDict

def extend_edge_atoms(extendEdgeAtomDict, originalBondList, validAtomSet):
    totalExtendedEdgeAtomList = []
    for atom, bondDist in extendEdgeAtomDict.items():
        # Shortcircuit this if bondDist is 0 - means the current edge atom is fine
        if bondDist == 0:
            totalExtendedEdgeAtomList.append(atom)

        # Find partial structure will work with a single edge atom
        extendedValidAtomSet, extendedEdgeAtomList = find_partial_structure(atom, originalBondList, None, bondDistance=bondDist)
        
        # Remove edge atoms that already existed in the validAtomSet - these will be atoms closer to the original bonding atoms
        extendedEdgeAtomList = [val for val in extendedEdgeAtomList if val not in validAtomSet]
        validAtomSet.update(extendedValidAtomSet)
        totalExtendedEdgeAtomList.extend(extendedEdgeAtomList)

    return validAtomSet, totalExtendedEdgeAtomList


def get_neighbours(atomIDList, bondsList, newBondAtoms):
    '''
    Get atomIDs of neighbouring atoms for each atom in atomIDList

    Bonding atoms don't appear as neighbours as this is used for symmetry checks.
    Bonding atoms always will have different neighbour fingerprints so no point looking at them.
    '''
    boundAtomsList = []

    # Determine what atoms are bound to an initial atom
    for atom in atomIDList:
        bondingAtoms = []
        for bond in bondsList:
            pairResult = pair_search(bond, atom)
            if pairResult is not None: # Stops bonding atoms appearing as neighbours CODE - and pairResult not in newBondAtoms
                bondingAtoms.append(pairResult)

        boundAtomsList.append([atom, bondingAtoms])

    # Create dictionary of initial atom keys and bound atom list values
    boundAtomsDict = {val[0]: val[1] for val in boundAtomsList}

    return boundAtomsDict

def get_additional_neighbours(neighboursDict, searchAtomID, searchNeighbours, bondingAtoms, unique=True):
    ''' Get atomIDs of the neighbours of a given atomID.     

        This is designed to get second and third neighbours of a given atomID. Further away
        neighbours are possible but may have unintended results.

        Args:
            unique: Prevent search from returning atomIDs that were already in the neighboursDict,
                in the searchNeighbours if specified, and the atomID. 

        Returns:
            List of neighbour atomIDs
    '''
   
    totalNeighbourSet = set()
    for currentNeighbour in searchNeighbours:
        totalNeighbourSet.update(neighboursDict[currentNeighbour])

    if unique:
        # Remove the original search atomID from totalNeighbourSet if present
        if searchAtomID in totalNeighbourSet:
            totalNeighbourSet.remove(searchAtomID)

        # Remove bonding atoms - don't want to use bonding atom fingerprints as they will always be different pre and post
        for bondingAtom in bondingAtoms:
            if bondingAtom in totalNeighbourSet:
                totalNeighbourSet.remove(bondingAtom)

        # Remove the neighbours from this search
        for currentNeighbour in searchNeighbours:
            if currentNeighbour in totalNeighbourSet:
                totalNeighbourSet.remove(currentNeighbour)
        
        # Remove initial neighbours from set if they aren't the searchNeighbours specified
        # This is for >= third neighbours
        if neighboursDict[searchAtomID] != searchNeighbours:
            for neighbour in neighboursDict[searchAtomID]:
                if neighbour in totalNeighbourSet:
                    totalNeighbourSet.remove(neighbour)

    return list(totalNeighbourSet)

def element_atomID_dict(fileName, elementsByType):
    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get charge
    data = clean_data(lines)
    sections = find_sections(data)
    try: # Try is for getting types from molecule file types
        types = get_data('Types', data, sections, useExcept=False)
    except ValueError: # Exception gets types from standard lammps file type
        atoms = get_data('Atoms', data, sections, useExcept=False)
        types = [[atomRow[0], atomRow[2]] for atomRow in atoms]
    typesDict = {row[0]: row[1] for row in types} # Keys: ID, Val: Type

    # Ensure elementsByType is uppercase
    elementsByTypeDict = {index+1: val.upper() for index, val in enumerate(elementsByType)} # Keys: Type, Val: Elements

    elementIDDict = {key: elementsByTypeDict[int(val)] for key, val in typesDict.items()}

    return elementIDDict
