import os
import math
import numpy as np
import matplotlib.pyplot as plt
from LammpsSearchFuncs import get_data, find_sections
from LammpsTreatmentFuncs import clean_data

def gen_distance_matrix(directory, fileName):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get coords
    data = clean_data(lines)
    sections = find_sections(data)
    coords = get_data('Coords', data, sections)

    # Collect atom IDs
    atomsIDs = [row[0] for row in coords]

    # Convert str coords to floats
    floatCoords = []
    for l in coords:
        floatRow = [float(val) for val in l]
        floatCoords.append(floatRow)

    # Calculate distance from atom1 to every other atom
    atom1Distances = []
    for atom1 in floatCoords:
        atom2Distances = []
        for atom2 in floatCoords:
            distance = math.sqrt((atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2 + (atom1[3] - atom2[3])**2)
            atom2Distances.append(distance)
        
        atom1Distances.append(atom2Distances)

    # Convert list of lists to matrix
    distMatrix = np.array(atom1Distances)

    return distMatrix, atomsIDs

def plot_matrix(matrixList):
    for matrix in matrixList:
        fig = plt.figure()
        ax = plt.gca()
        im = ax.matshow(matrix[0])
        fig.colorbar(im)
        ax.set_xticks(np.arange(len(matrix[0])))
        ax.set_xticklabels(matrix[1])
        ax.set_yticks(np.arange(len(matrix[0])))
        ax.set_yticklabels(matrix[1])
        
        
        # plt.matshow(matrix)
        # plt.colorbar()
      
    
    plt.show()

pre_matrix, pre_atomIDs = gen_distance_matrix('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction/', 'new_start_molecule.data')
post_matrix, post_atomIDs = gen_distance_matrix('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction/', 'new_post_rx1_molecule.data')


# Plot all matrices at the same time
# plot_matrix([[pre_matrix, pre_atomIDs], [post_matrix, post_atomIDs], [pre_matrix - post_matrix, pre_atomIDs]])

# diffMatrix = pre_matrix - post_matrix

searchRow = pre_matrix[11]
searchRowIndex = np.argsort(searchRow)
searchRowSorted = np.take_along_axis(searchRow, searchRowIndex, axis=0)
mask = [index for index, val in enumerate(searchRowSorted) if val > 10]
# for index in mask:
#     searchRowSorted[index] = 0

distArrayList = []

distDifference = []
indexList = []

for row in post_matrix:
    # Order row smallest to largest and store index
    rowIndex = np.argsort(row)
    indexList.append(rowIndex)

    rowSorted = np.take_along_axis(row, rowIndex, axis=0)
    # for index in mask:
    #     rowSorted[index] = 0

    # searchRow minus row
    diffRow = searchRowSorted - rowSorted
    distArrayList.append(np.around(diffRow, decimals=3))
    diffRow = np.abs(diffRow)
    diffRow = np.log(diffRow)
    diffRow = np.nan_to_num(diffRow, neginf=0.0)

    # Sum values
    arraySum = np.sum(diffRow)

    # Append
    distDifference.append(arraySum)

print()
# lowDiffIndex = []

# for searchRow in pre_matrix:
#     searchRowIndex = np.argsort(searchRow)
#     searchRowSorted = np.take_along_axis(searchRow, searchRowIndex, axis=0)
#     mask = [index for index, val in enumerate(searchRowSorted) if val > 10]
#     for index in mask:
#         searchRowSorted[index] = 0
    
    
#     distArrayList = []
    
#     distDifference = []
#     indexList = []

#     for row in post_matrix:
#         # Order row smallest to largest and store index
#         rowIndex = np.argsort(row)
#         indexList.append(rowIndex)

#         rowSorted = np.take_along_axis(row, rowIndex, axis=0)
#         for index in mask:
#             rowSorted[index] = 0

#         # searchRow minus row
#         diffRow = searchRowSorted - rowSorted
#         distArrayList.append(np.around(diffRow, decimals=2))
#         diffRow = np.abs(diffRow)

#         # Sum values
#         arraySum = np.sum(diffRow)

#         # Append
#         distDifference.append(arraySum)

#     # Find index of smallest value
#     _, smallestValIndex = min((val, idx) for (idx, val) in enumerate(distDifference))
    
#     lowDiffIndex.append(indexList[smallestValIndex])

# from statistics import mode
# finalResult = []
# for index in range(0,len(lowDiffIndex)-1):
#     indexValues = []
#     for l in lowDiffIndex:
#         indexValues.append(l[index])
    
#     modeVal = mode(indexValues)
#     finalResult.append(modeVal)

print()