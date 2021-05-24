import os
import math
import numpy as np
from natsort import natsorted
from collections import Counter, deque
from functools import reduce

from LammpsSearchFuncs import get_data, get_top_comments, read_top_comments, find_sections, pair_search, get_neighbours
from LammpsTreatmentFuncs import clean_data
from MappingFunctions import element_atomID_dict, calc_angles

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

class Atom():
    def __init__(self, atomID, element, bondingAtom, edgeAtom, NeighbourIDs):
        self.atomID = atomID
        self.element = element
        self.bondingAtom = bondingAtom
        self.edgeAtom = edgeAtom
        self.NeighbourIDs = NeighbourIDs
        self.NeighbourElements = []

    def check_mapped(self, mappedIDs, searchIndex):
        """Update NeighbourIDs.

        Updates NeighbourIDs by removing IDs that have already been mapped.
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
        
        self.NeighbourIDs = [ID for ID in self.NeighbourIDs if ID not in searchIndexMappedIDs]

    def get_elements(self, elementDict):
        """Gets elements for neighbour IDs.

        This is expected to be run after check_mapped, so that elements are only found
        for true  IDs.

        Args:
            elementDict: Dictionary of atomID keys, and element string values

        Returns:
            Updates existing class variable self.NeighbourElements
        """
        self.NeighbourElements = [elementDict[atomID]for atomID in self.NeighbourIDs]

    def map_elements(self, atomObject):
        """
        TO DO: How to handle multiple different atoms between base and target data?
            e.g. base and target are the same length but pre is H O H C and post is H C H C
        """

        # Output variables
        mapList = []
        missingPreAtoms = []
        queueAtoms = []

        # Match Function
        def matchNeighbour(preAtom, postAtom, preAtomIndex, postAtomIndex, mapList, queueList):
            mapList.append([preAtom.NeighbourIDs[preAtomIndex], postAtom.NeighbourIDs[postAtomIndex]])
            
            if self.NeighbourElements[preAtomIndex] != 'H':
                queueList.append([preAtom.NeighbourIDs[preAtomIndex], postAtom.NeighbourIDs[postAtomIndex]])

            postAtom.NeighbourIDs.pop(postAtomIndex)
            postAtom.NeighbourElements.pop(postAtomIndex)

        # Loop through neighbours for atom in one state and compare to neighbours of atom in other state
        for preIndex, neighbour in enumerate(self.NeighbourElements):
            elementOccurence = atomObject.NeighbourElements.count(neighbour)

            # If no match in post atom list it is a missingPreAtom
            if elementOccurence == 0:
                missingPreAtoms.append(self.NeighbourIDs[preIndex])
            
            # Assign atomIDs if there is only one matching element - could this go wrong if an element moves and an identical element takes its place?
            elif elementOccurence == 1:
                postIndex = atomObject.NeighbourElements.index(neighbour)
                matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)

            # More than one matching element requires additional considerations
            elif elementOccurence > 1:
                if neighbour == 'H': # H can be handled simply as all H are equivalent to each other in this case - ignores chirality
                    postHydrogenIndexList = [index for index, element in enumerate(atomObject.NeighbourElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    
                else:
                    # To come
                    # If looking at neighbours for the next atom, possible problem caused by popping and reducing ID and element lists
                    print("Symmetry atoms will be handled later")

        # Search mapList for missingPostAtoms
        mappedPostAtomList = [row[1] for row in mapList]
        missingPostAtoms = [neighbour for neighbour in atomObject.NeighbourIDs if neighbour not in mappedPostAtomList]

        return mapList, missingPreAtoms, missingPostAtoms, queueAtoms

def build_atom_objects(directory, fileName, elementDict):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get coords and bonds
    data = clean_data(lines)
    sections = find_sections(data)
    coords = get_data('Coords', data, sections)

    atomIDs = [row[0] for row in coords]
    bonds = get_data('Bonds', data, sections)

    # Get top comment info for bonding and edge atoms
    topComments = get_top_comments(lines)
    bondingAtoms = read_top_comments(topComments, 'Bonding_Atoms')
    edgeAtoms = read_top_comments(topComments, 'Edge_Atoms')
    edgeAtomFingerprints = read_top_comments(topComments, 'Edge_Atom_Fingerprints')
    
    # Check for if edge atoms are present
    if edgeAtoms is not None:
        edgeAtomDict = {edgeAtom: edgeAtomFingerprints[index] for index, edgeAtom in enumerate(edgeAtoms)}
    else:
        edgeAtomDict = None
        edgeAtoms = [] # Empty list passes later 'in' check

    # Build neighbours dict
    neighboursDict = get_neighbours(atomIDs, bonds)

    atomObjectList = []
    for atomID in atomIDs:
        neighbours = neighboursDict[atomID]

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

        atom = Atom(atomID, elementDict[atomID], bondingAtom, edgeAtom, neighbours)
        atomObjectList.append(atom)
    
    return atomObjectList, bondingAtoms, edgeAtomDict

# Returns atom class object that has specific atom ID
def get_atom_object(atomID, atomList):
    '''
    Returns the atom object for a specific atom ID. Potential uses include for turning queue
    atom IDs into atom objects for next stage of queue processing.
    Args:
        atomID: An integer atom ID
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


def map_from_path(directory, preFileName, postFileName, elementsByType):
    # Build atomID to element dict
    preElementDict = element_atomID_dict(directory, preFileName, elementsByType)
    postElementDict = element_atomID_dict(directory, postFileName, elementsByType)
    elementDictList = [preElementDict, postElementDict]

    # Generate atom class objects list
    preAtomObjectList, preBondingAtoms, preEdgeAtomDict = build_atom_objects(directory, preFileName, preElementDict)
    postAtomObjectList, postBondingAtoms, postEdgeAtomDict = build_atom_objects(directory, postFileName, postElementDict)

    # Initialise lists
    mappedIDList = []
    missingPreAtomsList = []
    missingPostAtomsList = []

    # Initialise queue
    queue = Queue()

    # Populate queue with bonding atoms and update mappedIDList
    queue_bond_atoms(preAtomObjectList, preBondingAtoms, postAtomObjectList, postBondingAtoms, mappedIDList, queue)

    # Populate queue with unique edge atoms if present
    queue_edge_atoms(preAtomObjectList, preEdgeAtomDict, postAtomObjectList, postEdgeAtomDict, mappedIDList, queue)

    # Search through queue creating new maps based on all elements in a given path
    while not queue.empty():
        currentAtoms = queue.get()
        for mainIndex, atom in enumerate(currentAtoms):
            atom.check_mapped(mappedIDList, mainIndex)
            atom.get_elements(elementDictList[mainIndex])
        
        newMap, missingPreAtoms, missingPostAtoms, queueAtoms = currentAtoms[0].map_elements(currentAtoms[1])

        # Convert queue atoms to atom class objects and add to queue
        add_to_queue(queue, queueAtoms, preAtomObjectList, postAtomObjectList)

        # Extend missing lists
        missingPreAtomsList.extend(missingPreAtoms)
        missingPostAtomsList.extend(missingPostAtoms)

        # Add new pairs to mapped ID list
        mappedIDList.extend(newMap)

    # Handle missing atoms
    # Update lists to remove any atoms that have been assigned and convert IDs into atom objects
    missingPreAtomsObjects = []
    missingPostAtomsObjects = []
    mappedPreAtoms = [pair[0] for pair in mappedIDList]
    mappedPostAtoms = [pair[1] for pair in mappedIDList]

    for atom in missingPreAtomsList:
        if atom not in mappedPreAtoms:
            atomsObject = get_atom_object(atom, preAtomObjectList)
            missingPreAtomsObjects.append(atomsObject)

    for atom in missingPostAtomsList:
        if atom not in mappedPostAtoms:
            atomsObject = get_atom_object(atom, postAtomObjectList)
            missingPostAtomsObjects.append(atomsObject)


    missingPostAtomElements = [atom.element for atom in missingPostAtomsObjects]
    for index, preAtom in enumerate(missingPreAtomsObjects):
        elementOccurence = missingPostAtomElements.count(preAtom.element)

        if elementOccurence == 0:
            print(f"Couldn't find a match in post missing atom for {preAtom.atomID}")

        elif elementOccurence == 1:
            postIndex = missingPostAtomElements.index(preAtom.element)
            mappedIDList.append([preAtom.atomID, missingPostAtomsObjects[postIndex].atomID])

        elif elementOccurence > 1:
            print(f"Too many element matches in post missing atom for {preAtom.atomID}")

    return mappedIDList
    
# map_from_path('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction', 'new_start_molecule.data', 'new_post_rx1_molecule.data', ['28', '62'], ['32', '15'], ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O'])


# Will take a lot more code in missing atoms to make it robust enough - possible revert points, multiple assignment of ID is possible currently

# Could add edge atom T/F to atom class, could help distinguish elements of same type

# Validation Checks:
# No repeated IDs / No IDs unassigned
# Ambiguous groups maintained pre and post aside from moved atoms
# Check that atom assigned from missingatoms is bound to the atoms that it expects to be bound to