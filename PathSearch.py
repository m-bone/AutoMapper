import os
import logging
from natsort import natsorted
from collections import Counter, deque

from LammpsSearchFuncs import get_data, get_top_comments, read_top_comments, find_sections, get_neighbours, get_additional_neighbours
from LammpsTreatmentFuncs import clean_data, save_text_file
from MappingFunctions import element_atomID_dict

# Classes and functions for search
class Queue:
    def __init__(self):
        self.elements = deque()
    
    def empty(self) -> bool:
        return not self.elements
    
    def add(self, x):
        for pair in x:
            self.elements.append(pair)
    
    def get(self):
        return self.elements.popleft()

def compare_symmetric_atoms(postNeighbourAtomObjectList, preNeighbourAtom, outputType):
    # Edge atom comparison

    # Neighbour comparison - no inference
    def compare_neighbours(neighbourLevel):
        neighbourComparison = [getattr(atomObject, neighbourLevel) for atomObject in postNeighbourAtomObjectList]
        neighbourFingerprint = [''.join(sorted(elements)) for elements in neighbourComparison] # sorted to get alphabetical fingerprints

        # Remove duplicate fingerprints
        countFingerprints = Counter(neighbourFingerprint)
        tuppledFingerprints = [(index, fingerprint) for index, fingerprint in enumerate(neighbourFingerprint) if countFingerprints[fingerprint] == 1]

        # Any of the potential post neighbours matches the pre atom fingerprint, return the post neighbour
        for index, fingerprint in tuppledFingerprints:
            if ''.join(getattr(preNeighbourAtom, neighbourLevel)) == fingerprint:
                logging.debug(f'Pre: {preNeighbourAtom.atomID}, Post: {postNeighbourAtomObjectList[index].atomID} found with {neighbourLevel}')
                if outputType == 'index':
                    return index
                elif outputType == 'atomID':
                    return postNeighbourAtomObjectList[index].atomID
                else:
                    print('Invalid output type specified for compare_symmetric_atoms')

    # First neighbour comparison
    symmetryResult = compare_neighbours('firstNeighbourElements')

    # Second neighbour comparison
    if symmetryResult is None:
        compare_neighbours('secondNeighbourElements')

    # Third neighbour comparison
    if symmetryResult is None:
        compare_neighbours('thirdNeighbourElements')

    # If it makes it through all these, guess assignment and warn user about this
    if symmetryResult is not None:
        return symmetryResult
    else:
        for index, postNeighbourAtom in enumerate(postNeighbourAtomObjectList):
            if postNeighbourAtom.element == preNeighbourAtom.element:
                logging.debug(f'Pre: {preNeighbourAtom.atomID}, Post: {postNeighbourAtom.atomID} found with symmetry inference')
                print(f'Note: Pre-bond atomID {preNeighbourAtom.atomID} has been assigned by inference to post-bond atomID {postNeighbourAtom.atomID}. Please check this is correct.')
                if outputType == 'index':
                        return index
                elif outputType == 'atomID':
                    return postNeighbourAtom.atomID
                else:
                    print('Invalid output type specified for compare_symmetric_atoms')


class Atom():
    def __init__(self, atomID, element, bondingAtom, edgeAtom, neighbourIDs, secondNeighbourIDs, thirdNeighbourIDs, neighbourElements, secondNeighbourElements, thirdNeighbourElements):
        self.atomID = atomID
        self.element = element
        self.bondingAtom = bondingAtom
        self.edgeAtom = edgeAtom

        # Neighbours
        self.mappedNeighbourIDs = neighbourIDs # This is changed according to mapping
        self.firstNeighbourIDs = neighbourIDs.copy() # This is fixed throughout mapping process
        self.secondNeighbourIDs = secondNeighbourIDs
        self.thirdNeighbourIDs = thirdNeighbourIDs

        self.mappedNeighbourElements = neighbourElements # This is changed according to mapping
        self.firstNeighbourElements = neighbourElements.copy() # This is fixed throughout mapping process
        self.secondNeighbourElements = secondNeighbourElements
        self.thirdNeighbourElements = thirdNeighbourElements

    def check_mapped(self, mappedIDs, searchIndex, elementDict):
        """Update neighbourIDs.

        Updates neighbourIDs by removing IDs that have already been mapped.
        This will be called before all neighbour mapping attempts to stop atoms
        being mapped multiple times.

        Args:
            mappedIDs: The total list of mappedIDs at this point in the mapping. This
                will contain pre- and post-atomIDs
            searchIndex: Determines whether to use pre- or post-atomIDs

        Returns:
            Updates existing class variable self.NeighbourIDs
        """
        searchIndexMappedIDs = [row[searchIndex] for row in mappedIDs]
        
        self.mappedNeighbourIDs = [ID for ID in self.mappedNeighbourIDs if ID not in searchIndexMappedIDs]
        self.mappedNeighbourElements = [elementDict[atomID]for atomID in self.mappedNeighbourIDs]


    def map_elements(self, atomObject, preAtomObjectList, postAtomObjectList):
        # Output variables
        mapList = []
        missingPreAtoms = []
        queueAtoms = []

        def allowed_maps(preAtom, postAtom):
            # Checks if elements appear the same number of times in pre and post atoms
            # If they don't, mapping is not allowed to take place and atoms are moved to missing lists
            preElementOccurences = Counter(preAtom.mappedNeighbourElements)
            postElementOccurences = Counter(postAtom.mappedNeighbourElements)

            allowedMapDict = {}
            for element, count in preElementOccurences.items():
                if count == postElementOccurences[element]:
                    allowedMapDict[element] = True
                else:
                    allowedMapDict[element] = False
                
            # Force all H to be be True as hydrogen can be mapped by inference
            if 'H' in allowedMapDict:
                allowedMapDict['H'] = True

            return allowedMapDict

        allowedMapDict = allowed_maps(self, atomObject)

        # Match Function
        def matchNeighbour(preAtom, postAtom, preAtomIndex, postAtomIndex, mapList, queueList):
            # Append pre and post atomIDs to map
            mapList.append([preAtom.mappedNeighbourIDs[preAtomIndex], postAtom.mappedNeighbourIDs[postAtomIndex]])
            
            # Add all non-hydrogen atom atomIDs to queue
            if preAtom.mappedNeighbourElements[preAtomIndex] != 'H':
                queueList.append([preAtom.mappedNeighbourIDs[preAtomIndex], postAtom.mappedNeighbourIDs[postAtomIndex]])

            # Remove post atomID from mappedID and mappedElement atom object values
            postAtom.mappedNeighbourIDs.pop(postAtomIndex)
            postAtom.mappedNeighbourElements.pop(postAtomIndex)

        # Loop through neighbours for atom in one state and compare to neighbours of atom in other state
        for preIndex, neighbour in enumerate(self.mappedNeighbourElements):
            elementOccurence = atomObject.mappedNeighbourElements.count(neighbour)

            # Check if maps with the neighbour element are allowed, if not add current element to missing list
            if allowedMapDict[neighbour] == False:
                missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])
                continue

            # If no match in post atom list it is a missingPreAtom
            if elementOccurence == 0:
                missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])
            
            # Assign atomIDs if there is only one matching element - could this go wrong if an element moves and an identical element takes its place?
            elif elementOccurence == 1:
                postIndex = atomObject.mappedNeighbourElements.index(neighbour)
                logging.debug(f'Pre: {self.mappedNeighbourIDs[preIndex]}, Post: {atomObject.mappedNeighbourIDs[postIndex]} found with single element occurence')
                matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)

            # More than one matching element requires additional considerations
            elif elementOccurence > 1:
                if neighbour == 'H': # H can be handled simply as all H are equivalent to each other in this case - ignores chirality
                    postHydrogenIndexList = [index for index, element in enumerate(atomObject.mappedNeighbourElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    logging.debug(f'Pre: {self.mappedNeighbourIDs[preIndex]}, Post: {atomObject.mappedNeighbourIDs[postIndex]} found with hydrogen symmetry inference')
                    matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    
                else:
                    # Get neighbour post atoms objects
                    postNeighbourIndices = [index for index, val in enumerate(atomObject.mappedNeighbourElements) if val == neighbour]
                    postNeighbourAtomIDs = [atomObject.mappedNeighbourIDs[i] for i in postNeighbourIndices]
                    postNeighbourAtomObjects = [postAtomObject for postAtomObject in postAtomObjectList if postAtomObject.atomID in postNeighbourAtomIDs]

                    # Get possible pre atom object
                    preNeighbourAtomObject = 0
                    for preAtomObject in preAtomObjectList:
                        if preAtomObject.atomID == self.mappedNeighbourIDs[preIndex]:
                            preNeighbourAtomObject = preAtomObject

                    # Find the post atom ID for the current pre atom
                    postNeighbourAtomID = compare_symmetric_atoms(postNeighbourAtomObjects, preNeighbourAtomObject, 'atomID')
                    if postNeighbourAtomID is not None:
                        postIndex = atomObject.mappedNeighbourIDs.index(postNeighbourAtomID)
                        matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    else:
                        # If no post atom found, add pre atom missing atom list
                        print(f'Could not find the symmetric pair for preAtom {self.mappedNeighbourIDs[preIndex]}')
                        missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])                                             


        # Search mapList for missingPostAtoms
        mappedPostAtomList = [row[1] for row in mapList]
        missingPostAtoms = [neighbour for neighbour in atomObject.mappedNeighbourIDs if neighbour not in mappedPostAtomList]

        return mapList, missingPreAtoms, missingPostAtoms, queueAtoms

def build_atom_objects(fileName, elementDict):
    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get coords and bonds
    data = clean_data(lines)
    sections = find_sections(data)
    coords = get_data('Coords', data, sections)

    atomIDs = [row[0] for row in coords]
    bonds = get_data('Bonds', data, sections)

    # Get top comment info for bonding, edge and delete atoms
    topComments = get_top_comments(lines)
    bondingAtoms = read_top_comments(topComments, 'Bonding_Atoms')
    edgeAtoms = read_top_comments(topComments, 'Edge_Atoms')
    edgeAtomFingerprints = read_top_comments(topComments, 'Edge_Atom_Fingerprints')
    deleteAtoms = read_top_comments(topComments, 'Delete_Atoms')
    
    # Check for if edge atoms are present
    if edgeAtoms is not None:
        edgeAtomDict = {edgeAtom: edgeAtomFingerprints[index] for index, edgeAtom in enumerate(edgeAtoms)}
    else:
        edgeAtomDict = None
        edgeAtoms = [] # Empty list passes later 'in' check

    # Build neighbours dict
    neighboursDict = get_neighbours(atomIDs, bonds)

    def get_elements(neighbourIDs, elementDict):
        return [elementDict[atomID]for atomID in neighbourIDs]

    atomObjectList = []
    for atomID in atomIDs:
        neighbours = neighboursDict[atomID]
        secondNeighbours = get_additional_neighbours(neighboursDict, atomID, neighbours)
        thirdNeighbours = get_additional_neighbours(neighboursDict, atomID, secondNeighbours)

        neighbourElements = get_elements(neighbours, elementDict)
        secondNeighbourElements = get_elements(secondNeighbours, elementDict)
        thirdNeighbourElements = get_elements(thirdNeighbours, elementDict)

        # Check if atom is a bonding atom, return boolean
        if atomID in bondingAtoms:
            bondingAtom = True
        else:
            bondingAtom = False

        # Check if atom is an edge atom, return fingerprint string or None
        if atomID in edgeAtoms:
            edgeAtomIndex = edgeAtoms.index(atomID)
            edgeAtom = edgeAtomFingerprints[edgeAtomIndex]
        else:
            edgeAtom = None

        atom = Atom(atomID, elementDict[atomID], bondingAtom, edgeAtom, neighbours, secondNeighbours, thirdNeighbours, neighbourElements, secondNeighbourElements, thirdNeighbourElements)
        atomObjectList.append(atom)
    
    return atomObjectList, bondingAtoms, edgeAtomDict, deleteAtoms

# Returns atom class object that has specific atom ID
def get_atom_object(atomID, atomList):
    '''
    Returns the atom object for a specific atom ID. Potential uses include turning queue
    atom IDs into atom objects for next stage of queue processing.
    Args:
        atomID: A string integer atom ID
        atomList: List of all possible atom objects, either pre- or post-bond

    Returns:
        Atom object
    '''
    for atom in atomList:
        if atom.atomID == atomID:
            return atom

def add_to_queue(queue, queueAtoms, preAtomObjectList, postAtomObjectList):
    queueAtomObjects = []
    for pair in queueAtoms:
        preAtom = get_atom_object(pair[0], preAtomObjectList)
        postAtom = get_atom_object(pair[1], postAtomObjectList)
        queueAtomObjects.append([preAtom, postAtom])
    queue.add(queueAtomObjects)

def queue_bond_atoms(preAtomObjectList, preBondingAtoms, postAtomObjectList, postBondingAtoms, mappedIDList, queue):
    # Loop through bonding atoms, getting atom objects and adding them to queue and mapped list
    for index, preBondAtom in enumerate(preBondingAtoms):
        preAtomObject = get_atom_object(preBondAtom, preAtomObjectList)
        postAtomObject = get_atom_object(postBondingAtoms[index], postAtomObjectList)
        queue.add([[preAtomObject, postAtomObject]])
        mappedIDList.append([preBondAtom, postBondingAtoms[index]])
        logging.debug(f'Pre: {preBondAtom}, Post: {postBondingAtoms[index]} found with user specified bond atom')


def queue_edge_atoms(preAtomObjectList, preEdgeAtomDict, postAtomObjectList, postEdgeAtomDict, mappedIDList, queue):
    # Skip function call if no edge atoms present
    if preEdgeAtomDict is None:
        return

    # Check edge atom fingerprints are the same
    assert natsorted(list(preEdgeAtomDict.values())) == natsorted(list(postEdgeAtomDict.values())), 'Pre and post edge atom fingerprints do not match.'
    
    # Find unique fingerprints 
    countEdgeAtomFingerprints = Counter(list(preEdgeAtomDict.values()))
    uniqueEdgeAtomFingerprints = [fingerprint for fingerprint in preEdgeAtomDict.values() if countEdgeAtomFingerprints[fingerprint] == 1]
    for fingerprint in uniqueEdgeAtomFingerprints:
        # Find index of fingerprint values, then get key at that index
        preEdgeAtomID = list(preEdgeAtomDict.keys())[list(preEdgeAtomDict.values()).index(fingerprint)]
        postEdgeAtomID = list(postEdgeAtomDict.keys())[list(postEdgeAtomDict.values()).index(fingerprint)]
        
        preAtomObject = get_atom_object(preEdgeAtomID, preAtomObjectList)
        postAtomObject = get_atom_object(postEdgeAtomID, postAtomObjectList)
        
        queue.add([[preAtomObject, postAtomObject]])
        mappedIDList.append([preEdgeAtomID, postEdgeAtomID])
        logging.debug(f'Pre: {preEdgeAtomID}, Post: {postEdgeAtomID} found with specified edge atom')

def map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList):
    # If delete atoms provided, add them to the mappedIDList. No purpose to including them in the queue
    if preDeleteAtoms is not None:
        assert postDeleteAtoms is not None, 'Delete atoms found in pre-bond file but not in post-bond file.'
        assert len(preDeleteAtoms) == len(postDeleteAtoms), 'Pre-bond and post-bond files have different numbers of delete atoms.'
        for index, preAtom in enumerate(preDeleteAtoms):
            mappedIDList.append([preAtom, postDeleteAtoms[index]])
            logging.debug(f'Pre: {preAtom}, Post: {postDeleteAtoms[index]} found with specified delete atom')

def run_queue(queue, mappedIDList, preAtomObjectList, postAtomObjectList, missingPreAtomList, missingPostAtomList, elementDictList):
    while not queue.empty():
        currentAtoms = queue.get()
        for mainIndex, atom in enumerate(currentAtoms):
            atom.check_mapped(mappedIDList, mainIndex, elementDictList[mainIndex])
        
        newMap, missingPreAtoms, missingPostAtoms, queueAtoms = currentAtoms[0].map_elements(currentAtoms[1], preAtomObjectList, postAtomObjectList)

        # Convert queue atoms to atom class objects and add to queue
        add_to_queue(queue, queueAtoms, preAtomObjectList, postAtomObjectList)

        # Extend missing lists
        missingPreAtomList.extend(missingPreAtoms)
        missingPostAtomList.extend(missingPostAtoms)

        # Add new pairs to mapped ID list
        mappedIDList.extend(newMap)

def get_missing_atom_objects(missingAtomList, atomObjectList):
    missingAtomObjects = []
    for atom in missingAtomList:
        atomObject = get_atom_object(atom, atomObjectList)
        missingAtomObjects.append(atomObject)

    return missingAtomObjects

def map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue):
    missingCheckCounter = 1
    while missingCheckCounter < 4 and len(missingPostAtomObjects) > 0:
        mappedPreAtomIndex = []
        for preIndex, preAtom in enumerate(missingPreAtomObjects):
            missingPostAtomElements = [atom.element for atom in missingPostAtomObjects]
            elementOccurence = missingPostAtomElements.count(preAtom.element)

            if elementOccurence == 0:
                print(f"Couldn't find a match in post missing atom for {preAtom.atomID}")

            elif elementOccurence == 1:
                postIndex = missingPostAtomElements.index(preAtom.element)
                logging.debug(f'Pre: {preAtom.atomID}, Post: {missingPostAtomObjects[postIndex].atomID} found with missing atoms single element occurence')
                matchMissing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                
            elif elementOccurence > 1:
                if preAtom.element == 'H':
                    postHydrogenIndexList = [index for index, element in enumerate(missingPostAtomElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    logging.debug(f'Pre: {preAtom.atomID}, Post: {missingPostAtomObjects[postIndex].atomID} found with missing atoms hydrogen symmetry inference')
                    matchMissing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                else:
                    potentialPostAtomObjects = [atomObject for atomObject in missingPostAtomObjects if atomObject.element == preAtom.element]
                    postIndex = compare_symmetric_atoms(potentialPostAtomObjects, preAtom, 'index')
                    if postIndex is not None:
                        matchMissing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                        logging.debug(f'The above atomID pair was found with missing atoms symmetry comparison')                                    

        # Refresh missingPreAtomObjects so that it doesn't print needless error messages on subsequent loops
        for index in sorted(mappedPreAtomIndex, reverse=True):
            del missingPreAtomObjects[index]        
        
        missingCheckCounter += 1

def matchMissing(preAtom, postAtomMissingIndex, missingPostAtomObjects, mapList, queue, preIndex, mappedPreAtomIndex):
    # Get post atom
    postAtom = missingPostAtomObjects[postAtomMissingIndex]

    mapList.append([preAtom.atomID, postAtom.atomID])
    
    if preAtom.element != 'H':
        queue.add([[preAtom, postAtom]]) # This circumvents add_to_queue()

    missingPostAtomObjects.pop(postAtomMissingIndex)
    mappedPreAtomIndex.append(preIndex)

def update_missing_list(missingAtomList, mappedIDList, mapIndex):
    mappedAtoms = [pair[mapIndex] for pair in mappedIDList]
    # Update missingAtomList to remove atoms that have been matched
    newMissingAtomList = [atom for atom in missingAtomList if atom not in mappedAtoms]

    return newMissingAtomList

def output_path(mappedIDList, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms):
    # Bonding atoms
    bondingAtoms = [['BondingIDs', '\n']]
    for atom in preBondingAtoms:
        bondingAtoms.extend([[atom]])
    bondingAtoms.extend(['\n'])

    # Delete Atoms
    deleteIDCount = []
    deleteAtoms = []
    if preDeleteAtoms is not None:
        deleteIDCount.extend([[str(len(preDeleteAtoms)) + ' deleteIDs']])
        deleteAtoms.extend([['DeleteIDs', '\n']])
        for atom in preDeleteAtoms:
            deleteAtoms.extend([[atom]])
        deleteAtoms.extend(['\n'])

    # Edge Atoms
    edgeIDCount = []
    edgeAtoms = []
    if preEdgeAtomDict is not None:
        edgeIDCount.extend([[str(len(preEdgeAtomDict)) + ' edgeIDs']])
        edgeAtoms.extend([['EdgeIDs', '\n']])
        for atom in preEdgeAtomDict.keys():
            edgeAtoms.extend([[atom]])
        edgeAtoms.extend(['\n'])
    edgeIDCount.extend('\n')

    # Equivalences
    equivalences = [['#This is an AutoMapper generated map\n'], [str(len(mappedIDList)) + ' equivalences']]
    equivalenceAtoms = [['Equivalences', '\n']]
    for atomPair in mappedIDList:
        equivalenceAtoms.extend([[atomPair[0] + '\t' + atomPair[1]]])

    # Output data
    output = []
    totalOutput = [equivalences, deleteIDCount, edgeIDCount, bondingAtoms, deleteAtoms, edgeAtoms, equivalenceAtoms]

    for section in totalOutput:
        output.extend(section)

    return output

def map_from_path(directory, preFileName, postFileName, elementsByType, logLevel):
    # Set log level
    logging.basicConfig(level=getattr(logging, logLevel.upper()))

    # Build atomID to element dict
    os.chdir(directory)
    preElementDict = element_atomID_dict(preFileName, elementsByType)
    postElementDict = element_atomID_dict(postFileName, elementsByType)
    elementDictList = [preElementDict, postElementDict]

    # Generate atom class objects list
    preAtomObjectList, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms = build_atom_objects(preFileName, preElementDict)
    postAtomObjectList, postBondingAtoms, postEdgeAtomDict, postDeleteAtoms = build_atom_objects(postFileName, postElementDict)

    # Initialise lists
    mappedIDList = []
    missingPreAtomList = []
    missingPostAtomList = []

    # Initialise queue
    queue = Queue()

    # Populate queue with bonding atoms and update mappedIDList
    queue_bond_atoms(preAtomObjectList, preBondingAtoms, postAtomObjectList, postBondingAtoms, mappedIDList, queue)

    # Populate queue with unique edge atoms if present
    queue_edge_atoms(preAtomObjectList, preEdgeAtomDict, postAtomObjectList, postEdgeAtomDict, mappedIDList, queue)

    # Add delete atoms to the mappedIDList
    map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList)

    # Search through queue creating new maps based on all elements in a given path
    run_queue(queue, mappedIDList, preAtomObjectList, postAtomObjectList, missingPreAtomList, missingPostAtomList, elementDictList)

    # Update missingPreAtoms to check if the missing atom search loop is needed
    missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)

    # If missing pre atoms are present, map missing atoms and rerun the queue until success or timeout
    timeoutCounter = 1
    while len(missingPreAtomList) > 0 and timeoutCounter < 6:
        # Update missing atom lists
        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        # Get pre atom objects
        missingPreAtomObjects = get_missing_atom_objects(missingPreAtomList, preAtomObjectList)

        # Add any post atoms that aren't in the map or already in the missing atoms
        mappedPostAtoms = [pair[1] for pair in mappedIDList]
        totalPostAtomList = [postAtom.atomID for postAtom in postAtomObjectList]
        unfoundMissingPostAtoms = [atomID for atomID in totalPostAtomList if atomID not in mappedPostAtoms and atomID not in missingPostAtomList]
        missingPostAtomList.extend(unfoundMissingPostAtoms)

        # Get post atom objects
        missingPostAtomObjects = get_missing_atom_objects(missingPostAtomList, postAtomObjectList)

        map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue)

        # Refresh missingAtomLists
        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        # Rerun the queue based on atom pairs added to queue from missingAtoms
        run_queue(queue, mappedIDList, preAtomObjectList, postAtomObjectList, missingPreAtomList, missingPostAtomList, elementDictList)
        print(f'missingPreAtoms after loop {timeoutCounter}: {missingPreAtomList}')

        timeoutCounter += 1

    # Order mappedIDList by preAtomID
    mappedIDList = natsorted(mappedIDList, key=lambda x: x[0])

    outputData = output_path(mappedIDList, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms)
    save_text_file('automap.data', outputData)

    return mappedIDList

# Validation Checks:
# No repeated IDs / No IDs unassigned
# Ambiguous groups maintained pre and post aside from moved atoms
# Check that atom assigned from missingatoms is bound to the atoms that it expects to be bound to