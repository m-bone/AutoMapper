import numpy as np
from sklearn.metrics import mean_absolute_error

from SpacialDistanceMatrix import spacial_distance_matrix
from BondDistanceMatrix import bond_distance_matrix
from ChargeMatrix import charge_matrix
from MappingFunctions import element_atomID_dict, element_validation

# File search constants
DATA_DIR = '/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction/'
PRE_FILE_NAME = 'new_start_molecule.data'
POST_FILE_NAME = 'new_post_rx1_molecule.data'
ELEMENT_BY_TYPE = ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O']
PRE_BONDING_ATOMS = ['28', '62']
POST_BONDING_ATOMS = ['32', '15']

# Get elements of atom IDs for pre and post molecules
preElements = element_atomID_dict(DATA_DIR, PRE_FILE_NAME, ELEMENT_BY_TYPE)
postElements = element_atomID_dict(DATA_DIR, POST_FILE_NAME, ELEMENT_BY_TYPE)

# Calculate bond distance matrices for pre and post molecule
preBondDistMat = bond_distance_matrix(DATA_DIR, PRE_FILE_NAME, PRE_BONDING_ATOMS)
postBondDistMat = bond_distance_matrix(DATA_DIR, POST_FILE_NAME, POST_BONDING_ATOMS) 

# Set value for hydrogen that moves to epoxide ring to zero - this will be automated / user supplied info in the future
for index, atomRow in enumerate(postBondDistMat):
    postBondDistMat[index][11] = 0.0

# Calculate spacial distance matrices for pre and post molecule
preDistMat, preAtomIDs = spacial_distance_matrix(DATA_DIR, PRE_FILE_NAME)
postDistMat, postAtomIDs = spacial_distance_matrix(DATA_DIR, POST_FILE_NAME)

# Get charge matrices for pre and post molecule
preChargeMat =  charge_matrix(DATA_DIR, PRE_FILE_NAME)
postChargeMat = charge_matrix(DATA_DIR, POST_FILE_NAME)

# Multiply bond distance and spacial distance matrices together
preBondMatrix = preBondDistMat   #* preDistMat * preChargeMat
postBondMatrix = postBondDistMat   #* postDistMat * postChargeMat

mappedIDList = []
for searchIndex, searchRow in enumerate(preBondMatrix):
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
        for row in postBondMatrix:
            # Sort row arrays from smallest to largest
            rowIndex = np.argsort(row)
            rowSorted = np.take_along_axis(row, rowIndex, axis=0)

            # MAE
            finalVal = mean_absolute_error(searchRowSorted, rowSorted)

            # Append - abs to get smallest value closest to zero
            distDifference.append(abs(finalVal))

        mappedPreAtomID, mappedPostAtomID = element_validation(preAtomIDs[searchIndex], postAtomIDs, distDifference, preElements, postElements)

        mappedIDList.append([mappedPreAtomID, mappedPostAtomID])
    
# Print test report
for mappedPair in mappedIDList:
    print(f'Atom {mappedPair[0]} is equivalent to atom {mappedPair[1]}')

correctPostAtomIDs = ['38', '39', '35', '41', '42', '32', '16', '5', '36', '37', '6', '4', '1', '3', '9', '17', '23', '15', '33', '34']
totalAtoms = len(correctPostAtomIDs)
correctAtoms = 0
incorrectPreAtomsList = []
for index, atom in enumerate(mappedIDList):
    if atom[1] == correctPostAtomIDs[index]:
        correctAtoms += 1
    else:
        incorrectPreAtomsList.append(atom[0])

print('\nTest Results:')
print(f'Correct atoms: {correctAtoms}. Accuracy: {round(correctAtoms / totalAtoms * 100, 1)}%')
print(f'Incorrect premolecule atomIDs: {incorrectPreAtomsList}')


# Next job is solve atom mapping for symmetrically equivalent atoms - i.e need 63 and 64 to have different mapped IDs
# 63 and 64 is an extra difficult case as 63 moves molecules - may require user input which will come in later versions