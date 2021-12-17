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
# A range of functions designed to search LAMMPS files for information.
# These functions work for 'read_data' files and 'molecule' files
##############################################################################
from natsort import natsorted
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
        
def edge_atom_fingerprint_ids(edgeAtomList, originalBondList, validAtomSet):
    # Get edge atom neighbours
    edgeAtomFingerprintDict = get_neighbours(edgeAtomList, originalBondList, []) # Bonding atoms given as blank list, edge atoms can never have bonding atoms as a neighbour so not a problem
    
    # Filter out validAtomIDs that are within the partial structure
    filteredFingerprintDict = {}
    for key, atomList in edgeAtomFingerprintDict.items():
        cutList = [atom for atom in atomList if atom not in validAtomSet]
        filteredFingerprintDict[key] = cutList

    return filteredFingerprintDict

def get_neighbours(atomIDList, bondsList):
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
            if pairResult is not None: # Stops bonding atoms appearing as neighbours
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

    # Assert that there are enough types in elementsByType for the highest type in the types variable
    largestType = int(natsorted(types, key=lambda x: x[1])[-1][1]) # Types are stored as lists of [AtomNumber, TypeNumber]
    assert len(elementsByType) >= largestType, 'EBT (elements by type) is missing values. Check that all types are present and separated with a space.'

    elementIDDict = {key: elementsByTypeDict[int(val)] for key, val in typesDict.items()}

    return elementIDDict

def get_header(tidiedData):
    '''
    Extract all the data from the header of a LAMMPS data file.
    Return a dictionary of keyword keys and listed numeric values
    '''
    
    # Find stop line by searching for first line starting with letters
    def get_stop_line():
        for index, line in enumerate(tidiedData):
            # Checks to get past the initial comment line(s):
            if index == 0: continue
            if line[0] == '#': continue

            if line[0].isalpha():
                return index

    headerStopLine = get_stop_line()

    # Build dictionary of header parts with keyword keys and list numeric values
    headerData = tidiedData[0:headerStopLine]
    headerDict = {'comment': []}
    for line in headerData:
        if line[0].isalpha() or line[0] == '#':
            headerDict['comment'].extend([line])
        else:
            # Break line by spaces
            cutLine = line.split()
            
            # Search through line to get the numeric values - list due to two box dimensions
            valueList = []
            keyList = []
            for element in cutLine:
                # Convert value to int, failing this a float, failing this skip it
                try:
                    valueList.append(int(element))
                except ValueError:
                    try:
                        valueList.append(float(element))
                    except ValueError:
                        keyList.append(element)

            # Create dict from assembled parts
            headerDict['_'.join(keyList)] = valueList
    
    return headerDict

def convert_header(header):
    '''Convert a header dictionary back to a list of lists of strings for output'''
    
    stringHeader = []
    for key, values in header.items():
        headerLine = [' '.join([str(val) for val in values])]
        if key != 'comment':
            headerLine.extend(key.split('_'))
        stringHeader.append(headerLine)
    
    return stringHeader