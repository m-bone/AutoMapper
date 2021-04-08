import os
import math
import numpy as np

from LammpsSearchFuncs import get_data, find_sections
from LammpsTreatmentFuncs import clean_data

def get_atomIDs(directory, fileName):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get coords
    data = clean_data(lines)
    sections = find_sections(data)
    coords = get_data('Coords', data, sections)

    # Collect atom IDs
    atomIDs = [row[0] for row in coords]

    return atomIDs

def element_atomID_dict(directory, fileName, elementsByType):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get charge
    data = clean_data(lines)
    sections = find_sections(data)
    try: # Try is for getting types from molecule file types
        types = get_data('Types', data, sections)
    except ValueError: # Exception gets types from standard lammps file type
        atoms = get_data('Atoms', data, sections)
        types = [[atomRow[0], atomRow[2]] for atomRow in atoms]
    typesDict = {row[0]: row[1] for row in types} # Keys: ID, Val: Type

    elementsByTypeDict = {index+1: val for index, val in enumerate(elementsByType)} # Keys: Type, Val: Elements

    elementIDDict = {key: elementsByTypeDict[int(val)] for key, val in typesDict.items()}

    return elementIDDict

def element_validation(preAtomID, postAtomIDList, differenceList, preElementDict, postElementDict, postBondingAtoms):
    # Make a copy of unchanged differenceList
    originalDifferenceList = differenceList.copy()

    # Find lowest difference post atom ID that is the same element as the pre atom ID
    checkElement = 1
    while checkElement:
        # Find smallest value and corresponding index
        val, idx = min((val, idx) for (idx, val) in enumerate(differenceList))
        # Find the smallest value's index in the original list
        originalIndex = originalDifferenceList.index(val)

        if postAtomIDList[originalIndex] in postBondingAtoms:
            # If chosen ID is one of the bondingAtoms, it's wrong so can be removed
            del differenceList[idx]
        elif preElementDict[preAtomID] == postElementDict[postAtomIDList[originalIndex]]:
            # If elements are the same return the pre and post atom IDs
            return preAtomID, postAtomIDList[originalIndex], originalIndex
        else:
            # If the elements are different delete the smallest value by index and try again
            del differenceList[idx]

def calc_angles(atomID, angleList, coordList):
    validAngles = [angle[2:] for angle in angleList if atomID in angle[2:]]
    coordDict = {row[0]: row[1:] for row in coordList}

    coordsForAngles = []
    for angle in validAngles:
        angleCoords = []
        for atom in angle:
            coordList = coordDict[atom]
            coordList = [float(coord) for coord in coordList]
            angleCoords.append(coordList)

        coordsForAngles.append(angleCoords)

    totalAngles = []
    for angle in coordsForAngles:
        a1 = np.array(angle[0])
        a2 = np.array(angle[1])
        a3 = np.array(angle[2])

        # Calculate magnitude of vectors - Calc bond length
        bondOne = math.sqrt((a1[0] - a2[0])**2 + (a1[1] - a2[1])**2 + (a1[2] - a2[2])**2)
        bondTwo = math.sqrt((a2[0] - a3[0])**2 + (a2[1] - a3[1])**2 + (a2[2] - a3[2])**2)

        # Adjust vectors for origin (the middle atom)
        vectorOne = a1 - a2
        vectorTwo = a3 - a2

        # Dot product of the vectors
        dotProd = np.dot(vectorOne, vectorTwo)

        # Angle calculation
        angle = dotProd / (bondOne*bondTwo)
        angle = np.degrees(np.arccos(angle))
        totalAngles.append(angle)
    
    return totalAngles