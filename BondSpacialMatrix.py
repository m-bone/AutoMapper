import numpy as np
from sklearn.metrics import mean_absolute_error
from SpacialDistanceMatrix import spacial_distance_matrix
from BondDistanceMatrix import bond_distance_matrix

# File search constants
DATA_DIR = '/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction/'
PRE_FILE_NAME = 'new_start_molecule.data'
POST_FILE_NAME = 'new_post_rx1_molecule.data'

# Calculate bond distance matrices for pre and post molecule
preBondDistMat = bond_distance_matrix(DATA_DIR, PRE_FILE_NAME, ['28', '62'])
postBondDistMat = bond_distance_matrix(DATA_DIR, POST_FILE_NAME, ['32', '15']) 

# Set value for hydrogen that moves to epoxide ring to zero - this will be automated / user supplied info in the future
for index, atomRow in enumerate(postBondDistMat):
    postBondDistMat[index][11] = 0.0

# Calculate spacial distance matrices fro pre and post molecule
preDistMat, preAtomIDs = spacial_distance_matrix(DATA_DIR, PRE_FILE_NAME)
postDistMat, postAtomIDs = spacial_distance_matrix(DATA_DIR, POST_FILE_NAME)

# Multiply bond distance and spacial distance matrices together
preBondMatrix = preBondDistMat * preDistMat
postBondMatrix = postBondDistMat * postDistMat

for searchIndex, searchRow in enumerate(preBondMatrix):
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

    # Determine index of smallest difference value
    val, idx = min((val, idx) for (idx, val) in enumerate(distDifference))

    # Current result output
    print(f'Atom {preAtomIDs[searchIndex]} is equivalent to atom {postAtomIDs[idx]}')