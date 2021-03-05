import os
import numpy as np
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
    for index, atomRow in enumerate(postBondDistMat):
        postBondDistMat[index][atomIndex] = 0.0

# Sample weight matrix - lower weight for atoms with significant movement
sampleWeights = np.ones(len(postAtomIDs))
for atom in POST_MINOR_MOVED_ATOMS:
    atomIndex = postAtomIDs.index(atom)
    sampleWeights[atomIndex] = _WEIGHT_COEFF

mappedIDList = []
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

        mappedPreAtomID, mappedPostAtomID = element_validation(preAtomIDs[searchIndex], postAtomIDs, distDifference, preElements, postElements)

        mappedIDList.append([mappedPreAtomID, mappedPostAtomID])



# Print test report
for mappedPair in mappedIDList:
    print(f'Atom {mappedPair[0]} is mapped to atom {mappedPair[1]}')

correctPostAtomIDs = [['38'], ['39'], ['35'], ['41', '42'], ['42', '41'], ['32'], ['16'], ['5', '36'], ['36', '5'], ['37'], ['6', '9'], ['4'], ['1', '3'], ['3', '1'], ['9', '6'], ['17', '23'], ['23', '17'], ['15'], ['33', '34'], ['34', '33']]
totalAtoms = len(correctPostAtomIDs)
correctAtoms = 0
incorrectPreAtomsList = []
for index, atom in enumerate(mappedIDList):
    if atom[1] in correctPostAtomIDs[index]:
        correctAtoms += 1
    else:
        incorrectPreAtomsList.append(atom[0])

print(f'Test Results: Weight coeff is {_WEIGHT_COEFF}')
print(f'Correct atoms: {correctAtoms}. Accuracy: {round(correctAtoms / totalAtoms * 100, 1)}%')
print(f'Incorrect premolecule atomIDs: {incorrectPreAtomsList}')