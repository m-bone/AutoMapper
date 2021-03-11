import os
import math
import numpy as np
from natsort import natsorted
from sklearn.metrics import mean_absolute_error

from BondDistanceMatrix import bond_distance_matrix
from MappingFunctions import element_atomID_dict, element_validation, get_atomIDs

# File search constants and user inputs
DATA_DIR = os.getcwd() + '/Test_Cases'
PRE_FILE_NAME = 'new_start_molecule.data'
POST_FILE_NAME = 'new_post_rx1_molecule.data'
ELEMENT_BY_TYPE = ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O']
PRE_BONDING_ATOMS = ['28', '62']
POST_BONDING_ATOMS = ['32', '15']
POST_MAJOR_MOVED_ATOMS = ['33']
POST_MINOR_MOVED_ATOMS = ['16']


def atom_mapping(DATA_DIR, PRE_FILE_NAME, POST_FILE_NAME, ELEMENT_BY_TYPE, PRE_BONDING_ATOMS, POST_BONDING_ATOMS, POST_MAJOR_MOVED_ATOMS, POST_MINOR_MOVED_ATOMS):
    _WEIGHT_COEFF = 0.0
    # Get elements of atom IDs for pre and post molecules
    preElements = element_atomID_dict(DATA_DIR, PRE_FILE_NAME, ELEMENT_BY_TYPE)
    postElements = element_atomID_dict(DATA_DIR, POST_FILE_NAME, ELEMENT_BY_TYPE)

    # Get atomIDs - using existing function for now 
    preAtomIDs = get_atomIDs(DATA_DIR, PRE_FILE_NAME)
    postAtomIDs = get_atomIDs(DATA_DIR, POST_FILE_NAME)

    # Calculate bond distance matrices for pre and post molecule
    preBondDistMat = bond_distance_matrix(DATA_DIR, PRE_FILE_NAME, PRE_BONDING_ATOMS)
    postBondDistMat = bond_distance_matrix(DATA_DIR, POST_FILE_NAME, POST_BONDING_ATOMS) 

    # Set value for hydrogen that moves to epoxide ring to zero - this will be automated / user supplied info in the future
    for atom in POST_MAJOR_MOVED_ATOMS:
        atomIndex = postAtomIDs.index(atom)
        for index, _ in enumerate(postBondDistMat):
            postBondDistMat[index][atomIndex] = 0.0

    # Sample weight matrix - lower weight for atoms with significant movement
    sampleWeights = np.ones(len(postAtomIDs))
    for atom in POST_MINOR_MOVED_ATOMS:
        atomIndex = postAtomIDs.index(atom)
        sampleWeights[atomIndex] = _WEIGHT_COEFF

    mappedIDList = []
    mappedPostAtomsIndex = []
    for searchIndex, searchRow in enumerate(preBondDistMat):
        # Shortcircuit this search if search atom is a bonding atom
        if preAtomIDs[searchIndex] in PRE_BONDING_ATOMS:
            bondingIndex = PRE_BONDING_ATOMS.index(preAtomIDs[searchIndex])
            bondingPostAtomID = POST_BONDING_ATOMS[bondingIndex]
            mappedIDList.append([preAtomIDs[searchIndex], bondingPostAtomID])

        else:
            # Sort search row arrays from smallest to largest
            searchRowIndex = np.argsort(searchRow)
            searchRowSorted = np.take_along_axis(searchRow, searchRowIndex, axis=0)
            
            distDifference = []
            for row in postBondDistMat:
                # Sort row arrays from smallest to largest
                rowIndex = np.argsort(row)
                rowSorted = np.take_along_axis(row, rowIndex, axis=0)
                
                # Sort sample weight matrix the same as row
                sampleWeightsSorted = np.take_along_axis(sampleWeights, rowIndex, axis=0)

                # MAE
                finalVal = mean_absolute_error(searchRowSorted, rowSorted, sample_weight=sampleWeightsSorted)

                # Append - abs to get smallest value closest to zero
                distDifference.append(abs(finalVal))

            mappedPreAtomID, mappedPostAtomID, postAtomIDIndex = element_validation(preAtomIDs[searchIndex], postAtomIDs, distDifference, preElements, postElements, POST_BONDING_ATOMS)

            mappedIDList.append([mappedPreAtomID, mappedPostAtomID])
            mappedPostAtomsIndex.append(postAtomIDIndex)

    # Ambiguous Atom Group Processing
    # Gather all the pairs
    mappedPostAtomIDs = [val[1] for val in mappedIDList]
    repeatedPostIDSet = natsorted(set([val for val in mappedPostAtomIDs if mappedPostAtomIDs.count(val) > 1]))
    repeatedIndexes = [postAtomIDs.index(ID) for ID in repeatedPostIDSet]

    ambiguousGroupPairs = []
    # Loop through all post atoms to find similar
    for index in repeatedIndexes:
        matchArray = postBondDistMat[index]
        
        # Sort search row arrays from smallest to largest
        searchRowIndex = np.argsort(matchArray)
        searchRowSorted = np.take_along_axis(matchArray, searchRowIndex, axis=0)
            
        distDifference = []
        for row in postBondDistMat:
            # Sort row arrays from smallest to largest
            rowIndex = np.argsort(row)
            rowSorted = np.take_along_axis(row, rowIndex, axis=0)
            
            # Sort sample weight matrix the same as row
            sampleWeightsSorted = np.take_along_axis(sampleWeights, rowIndex, axis=0)

            # MAE
            finalVal = mean_absolute_error(searchRowSorted, rowSorted, sample_weight=sampleWeightsSorted)

            # Append - abs to get smallest value closest to zero
            distDifference.append(abs(finalVal))

        # Set repeatedIndex value to nan as it will always be zero
        distDifference[index] = math.nan
        _, smallestIndex = min((val, idx) for (idx, val) in enumerate(distDifference))

        ambiguousGroupPairs.append([postAtomIDs[index], postAtomIDs[smallestIndex]])
        # print(f'Atom {postAtomIDs[index]} is paired to atom {postAtomIDs[smallestIndex]}')

    # Update mappedIDList based on the ambiguousGroupPairs values
    # Interestingly, mappedIDList can be updated with the iterator, but ambiguousGroupPairs needs to be deleted with the index value
    for mappedID in mappedIDList:
        if mappedID[1] in repeatedPostIDSet: # If mappedPostAtomID is one that is repeated
            for index, groupPair in enumerate(ambiguousGroupPairs):
                if groupPair[0] == mappedID[1]: # If groupPair is a matching PostAtomID
                    mappedID[1] = groupPair[1]
                    del ambiguousGroupPairs[index]
                    break

    return mappedIDList
# This needs to include a check to make sure it's not updating the atom to an ID already assigned - might be best to have an unassigned-postAtomIDList

