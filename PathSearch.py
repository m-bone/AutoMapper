import os
import math
import numpy as np
from natsort import natsorted
from collections import deque
from functools import reduce

from LammpsSearchFuncs import get_data, find_sections, pair_search
from LammpsTreatmentFuncs import clean_data
from MappingFunctions import element_atomID_dict, calc_angles

# Classes and functions for search
class Graph:
    def __init__(self):
        self.edges = {}
    
    def neighbours(self, id):
        return self.edges[id]

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

def fill_graph(atomIDList, bondsList):
    # Load graph object
    moleculeGraph = Graph()
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
    moleculeGraph.edges = boundAtomsDict

    return moleculeGraph

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

class Atom():
    def __init__(self, atomID, bondingAtom, NeighbourIDs):
        self.atomID = atomID
        self.bondingAtom = bondingAtom
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
        def match(preAtom, postAtom, preAtomIndex, postAtomIndex, mapList, queueList):
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
                match(self, atomObject, preIndex, postIndex, mapList, queueAtoms)

            # More than one matching element requires additional considerations
            elif elementOccurence > 1:
                if neighbour == 'H': # H can be handled simply as all H are equivalent to each other in this case - ignores chirality
                    postHydrogenIndexList = [index for index, element in enumerate(atomObject.NeighbourElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    match(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    
                else:
                    # To come
                    # If looking at neighbours for the next atom, possible problem caused by popping and reducing ID and element lists
                    print("Symmetry atoms will be handled later")

        # Search mapList for missingPostAtoms
        mappedPostAtomList = [row[1] for row in mapList]
        missingPostAtoms = [neighbour for neighbour in atomObject.NeighbourIDs if neighbour not in mappedPostAtomList]

        return mapList, missingPreAtoms, missingPostAtoms, queueAtoms

def build_atom_objects(directory, fileName, bondingAtoms, elementsByType):
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

    # Build neighbours dict
    neighboursDict = get_neighbours(atomIDs, bonds)

    atomList = []
    for atom in atomIDs:
        neighbours = neighboursDict[atom]
        if atom in bondingAtoms:
            bondingAtom = True
        else:
            bondingAtom = False
        atom = Atom(atom, bondingAtom, neighbours)
        atomList.append(atom)
    
    return atomList

# Returns atom class object that has specific atom ID
def get_atom_object(atomID, atomList):
    for atom in atomList:
        if atom.atomID == atomID:
            return atom

def map_from_path(directory, preFileName, postFileName, preBondingAtoms, postBondingAtoms, elementsByType):
    preAtomObjectList = build_atom_objects(directory, preFileName, preBondingAtoms, elementsByType)
    postAtomObjectList = build_atom_objects(directory, postFileName, postBondingAtoms, elementsByType)

    # Build atomID to element dict
    preElementDict = element_atomID_dict(directory, preFileName, elementsByType)
    postElementDict = element_atomID_dict(directory, postFileName, elementsByType)
    elementDictList = [preElementDict, postElementDict]

    # Initialise lists
    mappedIDList = []
    missingPreAtomsList = []
    missingPostAtomsList = []

    # Add bonding atoms to mapped ID list
    for index, atom in enumerate(preBondingAtoms):
        mappedIDList.append([atom, postBondingAtoms[index]])

    # Mapping loop
    for bondingIndex, preBondAtom in enumerate(preBondingAtoms):
        preStartAtom = get_atom_object(preBondAtom, preAtomObjectList)
        postStartAtom = get_atom_object(postBondingAtoms[bondingIndex], postAtomObjectList)
        
        queue = Queue()
        queue.add([[preStartAtom, postStartAtom]])

        # Further iterations of this loop are going to need to queue with atomIDs so the correct pre and post neighbours can be compared
        # Would fail if it tried to compare two atom centres that aren't the confirmed equivalent map

        while not queue.empty():
            currentAtoms = queue.get()
            for mainIndex, atom in enumerate(currentAtoms):
                atom.check_mapped(mappedIDList, mainIndex)
                atom.get_elements(elementDictList[mainIndex])
            
            newMap, missingPreAtoms, missingPostAtoms, queueAtoms = currentAtoms[0].map_elements(currentAtoms[1])

            # Convert queue atoms to atom class objects and to queue
            queueAtomObjects = []
            for pair in queueAtoms:
                preAtom = get_atom_object(pair[0], preAtomObjectList)
                postAtom = get_atom_object(pair[1], postAtomObjectList)
                queueAtomObjects.append([preAtom, postAtom])
            queue.add(queueAtomObjects)

            # Extend missing lists
            missingPreAtomsList.extend(missingPreAtoms)
            missingPostAtomsList.extend(missingPostAtoms)

            # Add new pairs to mapped ID list
            mappedIDList.extend(newMap)


    print()


map_from_path('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction', 'new_start_molecule.data', 'new_post_rx1_molecule.data', ['28', '62'], ['32', '15'], ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O'])

# Could add edge atom T/F to atom class, could help distinguish elements of same type

# Validation Checks:
# No repeated IDs / No IDs unassigned
# Ambiguous groups maintained pre and post aside from moved atoms