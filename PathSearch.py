import os
import math
import numpy as np
from natsort import natsorted
from collections import deque
from functools import reduce

from LammpsSearchFuncs import get_data, find_sections, pair_search
from LammpsTreatmentFuncs import clean_data
from MappingFunctions import element_atomID_dict

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
    
    def put(self, x):
        self.elements.append(x)
    
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
        self.NeighbourElements = [elementDict[atomID] for atomID in self.NeighbourIDs]

    def map_elements(self, atomObject):
        # Output variables
        mapList = []
        missingPreAtoms = []
        missingPostAtoms = []

        # Hydrogen handling variables
        postHydrogenIDList = []
        createHydrogenList = True

        for index, neighbour in enumerate(self.NeighbourElements):
            # Count from the longest neighbour list, but use the atomObject list if both lists are equal length
            if len(self.NeighbourElements) > len(atomObject.NeighbourElements):              
                elementOccurence = self.NeighbourElements.count(neighbour)
                # primaryData = atomObject
                # secondaryData = self
            else:
                elementOccurence = atomObject.NeighbourElements.count(neighbour)
                # primaryData = self
                # secondaryData = atomObject

            if elementOccurence == 0:
                continue
            elif elementOccurence == 1:
                matchIndex = atomObject.NeighbourElements.index(neighbour)
                mapList.append([self.NeighbourIDs[index], atomObject.NeighbourIDs[matchIndex]])
            elif elementOccurence > 1:
                if neighbour == 'H':
                    if createHydrogenList:
                        createHydrogenList = False
                        postHydrogenIDList = [atomObject.NeighbourIDs[index] for index, element in enumerate(atomObject.NeighbourElements) if element == 'H']
                    elif createHydrogenList == False and len(postHydrogenIDList) == 0:
                        missingPreAtoms.append(self.NeighbourIDs[index])
                        continue

                    mapList.append([self.NeighbourIDs[index], postHydrogenIDList.pop()])
                else:
                    print("AHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH")

        return mapList, missingPreAtoms

def path_search(directory, fileName, bondingAtoms, elementsByType):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get coords and bonds
    data = clean_data(lines)
    sections = find_sections(data)
    coords = get_data('Coords', data, sections)

    coordDict = {row[0]: row[1:] for row in coords}

    atomIDs = [row[0] for row in coords]
    bonds = get_data('Bonds', data, sections)

    # Build neighbours dict
    neighboursDict = get_neighbours(atomIDs, bonds)

    atomList = []
    for bondAtom in bondingAtoms:
        neighbours = neighboursDict[bondAtom]
        
        atom = Atom(bondAtom, True, neighbours)
        atomList.append(atom)
    
    return atomList

def find_start_atom(atomList, startAtom):
    for atom in atomList:
        if atom.atomID == startAtom:
            return atom

def map_from_path(directory, preFileName, postFileName, preBondingAtoms, postBondingAtoms, elementsByType):
    preBond = path_search(directory, preFileName, preBondingAtoms, elementsByType)
    postBond = path_search(directory, postFileName, postBondingAtoms, elementsByType)

    # Build atomID to element dict
    preElementDict = element_atomID_dict(directory, preFileName, elementsByType)
    postElementDict = element_atomID_dict(directory, postFileName, elementsByType)
    elementDictList = [preElementDict, postElementDict]

    mappedIDList = []

    # Add bonding atoms to mapped ID list
    for index, atom in enumerate(preBondingAtoms):
        mappedIDList.append([atom, postBondingAtoms[index]])

    # Mapping loop
    for bondingIndex, preBondAtom in enumerate(preBondingAtoms):
        preStartAtom = find_start_atom(preBond, preBondAtom)
        postStartAtom = find_start_atom(postBond, postBondingAtoms[bondingIndex])
        
        queue = Queue()
        queue.put([preStartAtom, postStartAtom])

        while not queue.empty():
            currentAtoms = queue.get()
            for mainIndex, atom in enumerate(currentAtoms):
                atom.check_mapped(mappedIDList, mainIndex)
                atom.get_elements(elementDictList[mainIndex])
            
            newMap, missingPreAtoms = currentAtoms[0].map_elements(currentAtoms[1])
            print()

def breadth_first_search(graph, start, target):    
    # This code and associated classes is adapted from from https://www.redblobgames.com/pathfinding/a-star/introduction.html
    # and https://www.redblobgames.com/pathfinding/a-star/implementation.html
    queue = Queue()
    queue.put(start)
    came_from = dict()
    came_from[start] = None
    
    while not queue.empty(): # If queue is not empty then there are atoms that haven't been searched
        current = queue.get() # Removes one from queue and looks at it
        for next in graph.neighbours(current): # Get all the neighbours of the atom
            if next not in came_from: # Neighbours haven't been looked at
                queue.put(next) # Add them to the queue
                came_from[next] = current # Add to dict

    current = target
    path = []
    while current != start:
        path.append(current)
        try:
            current = came_from[current]
        except KeyError: # If no path exists, path is None
            path = [None]
            break
    if path[0] is not None:
        path.append(start)
        path.reverse()

    return path

map_from_path('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction', 'new_start_molecule.data', 'new_post_rx1_molecule.data', ['28', '62'], ['32', '15'], ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O'])