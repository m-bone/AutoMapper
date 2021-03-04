import os
import math
import numpy as np
from natsort import natsorted
from collections import deque
from functools import reduce

from LammpsSearchFuncs import get_data, find_sections, pair_search
from LammpsTreatmentFuncs import clean_data

def calc_bond_length(bondRow, coordDict):
    # Get IDs
    bondID = bondRow[0]
    atom1ID = bondRow[2]
    atom2ID = bondRow[3]

    # Get coords
    atom1Coords = coordDict[atom1ID]
    atom2Coords = coordDict[atom2ID]

    # Convert coords to floats
    atom1Coords = [float(coord) for coord in atom1Coords]
    atom2Coords = [float(coord) for coord in atom2Coords]

    # Calculate bond length
    bondLength = math.sqrt((atom1Coords[0] - atom2Coords[0])**2 + (atom1Coords[1] - atom2Coords[1])**2 + (atom1Coords[2] - atom2Coords[2])**2)

    return [bondID, bondLength]

def get_bond_path(atomList, bonds):
    bondIDList = []

    # If no bond is present, return empty list
    if atomList[0] is None:
        return bondIDList
    
    for index, atom in enumerate(atomList):
        # Get bond partner or quit if no more partners
        try:
            bondPartner = atomList[index+1]
        except IndexError:
            break

        # Sort atoms smallest to largest, this is a LAMMPS assumption
        bondPair = natsorted([atom, bondPartner])

        # Find bond that features both atoms
        for bond in bonds:
            if bond[2] == bondPair[0] and bond[3] == bondPair[1]:
                bondIDList.append(bond[0])

    return bondIDList

def calc_path_distance(bondList, bondDict):
    # If bondList is empty return zero
    if len(bondList) == 0:
        return 0.0

    bondDistList = []
    for bondID in bondList:
        bondDistList.append(bondDict[bondID])

    bondDistMultiple = reduce((lambda x, y: x * y), bondDistList)
    return bondDistMultiple

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

def breadth_first_search(graph, start, target):
    # Shortcircuit to None path list if start and target are the same atom
    if start == target:
        return [None]
    
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

def bond_distance_matrix(directory, fileName, bondingAtoms):
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

    # Calculate bond lengths and build in a dictionary
    bondsLengthList = [calc_bond_length(bond, coordDict) for bond in bonds]
    bondLengthDict = {bond[0]: bond[1] for bond in bondsLengthList}

    # Edit bondLengthDict to make bonding atoms bond length zero
    sortBondingAtoms = natsorted(bondingAtoms)
    bondingAtomsBondID = []
    for bond in bonds:
        if bond[2] == sortBondingAtoms[0] and bond[3] == sortBondingAtoms[1]:
            bondingAtomsBondID.append(bond[0])
    if len(bondingAtomsBondID) > 0: # Stops this happening in pre-molecule cases before bond occurs
        bondLengthDict[bondingAtomsBondID[0]] = 0.0

    # Create graph of bonding atom pairs for bond path search
    moleculeGraph = fill_graph(atomIDs, bonds) 

    # For a starting atom calculate the path to all atoms, get the bond IDs that make this path
    # calculate the bond distances and append to form a list of lists
    totalBondDistanceList = []
    for startAtom in atomIDs:
        atomBondDistanceList = []
        for otherAtom in atomIDs:
            atomPath = breadth_first_search(moleculeGraph, startAtom, otherAtom)
            bondPath = get_bond_path(atomPath, bonds)
            pathDistance = calc_path_distance(bondPath, bondLengthDict)
            atomBondDistanceList.append(pathDistance)
        
        totalBondDistanceList.append(atomBondDistanceList)

    # Convert list of lists to numpy matrix
    bondDistanceMatrix = np.array(totalBondDistanceList)
    
    return bondDistanceMatrix



