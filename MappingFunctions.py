import os
import numpy as np

from LammpsSearchFuncs import get_data, find_sections
from LammpsTreatmentFuncs import clean_data

def element_atomID_dict(directory, fileName, elementsByType):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get charge
    data = clean_data(lines)
    sections = find_sections(data)
    types = get_data('Types', data, sections)
    typesDict = {row[0]: row[1] for row in types} # Keys: ID, Val: Type

    elementsByTypeDict = {index+1: val for index, val in enumerate(elementsByType)} # Keys: Type, Val: Elements

    elementIDDict = {key: elementsByTypeDict[int(val)] for key, val in typesDict.items()}

    return elementIDDict

def element_validation(preAtomID, postAtomIDList, differenceList, preElementDict, postElementDict):
    # Make a copy of unchanged differenceList
    originalDifferenceList = differenceList.copy()

    checkElement = 1

    # Find lowest difference post atom ID that is the same element as the pre atom ID
    while checkElement:
        # Find smallest value and corresponding index
        val, idx = min((val, idx) for (idx, val) in enumerate(differenceList))
        # Find the smallest value's index in the original list
        originalIndex = originalDifferenceList.index(val)

        # If elements are the same return the pre and post atom IDs 
        if preElementDict[preAtomID] == postElementDict[postAtomIDList[originalIndex]]:
            return preAtomID, postAtomIDList[originalIndex]
        else:
            # If the elements are different delete the smallest value by index and try again
            del differenceList[idx]
