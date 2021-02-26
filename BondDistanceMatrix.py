import os
import math
import numpy as np
from collections import deque
import matplotlib.pyplot as plt
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

# def calc_bond_distance(startAtom, targetAtom, bonds, bondLengthDict):
#     searchPoint = startAtom
#     possibleStartList = []
#     for bond in bonds:
#         pairResult = pair_search(bond, searchPoint)
#         if pairResult is not None:
#             possibleStartList.append(pairResult)

class Graph:
    def __init__(self):
        self.edges = {}
    
    def neighbours(self, id):
        return self.edges[id]


def fill_graph(atomIDList, bondsList):
    moleculeGraph = Graph()
    boundAtomsList = []

    for atom in atomIDList:
        bondingAtoms = []
        for bond in bondsList:
            pairResult = pair_search(bond, atom)
            if pairResult is not None:
                bondingAtoms.append(pairResult)

        boundAtomsList.append([atom, bondingAtoms])

    boundAtomsDict = {val[0]: val[1] for val in boundAtomsList}
    moleculeGraph.edges = boundAtomsDict

    return moleculeGraph

class Queue:
    def __init__(self):
        self.elements = deque()
    
    def empty(self) -> bool:
        return not self.elements
    
    def put(self, x):
        self.elements.append(x)
    
    def get(self):
        return self.elements.popleft()

def breadth_first_search(graph, start, target):
    # This code and associated classes comes from https://www.redblobgames.com/pathfinding/a-star/introduction.html
    # and https://www.redblobgames.com/pathfinding/a-star/implementation.html
    queue = Queue()
    queue.put(start)
    came_from = dict()
    came_from[start] = None
    
    while not queue.empty():
        current = queue.get()
        for next in graph.neighbours(current):
            if next not in came_from:
                queue.put(next)
                came_from[next] = current

    current = target
    path = []
    while current != start:
        path.append(current)
        current = came_from[current]
    path.append(start)
    path.reverse()

    print(path)

def bond_distance_matrix(directory, fileName):
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

    moleculeGraph = fill_graph(atomIDs, bonds) 


    print('Reachable from: 22')
    breadth_first_search(moleculeGraph, '22', '34')

    print()

bond_distance_matrix('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction/', 'new_start_molecule.data')