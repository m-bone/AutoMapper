##############################################################################
# Developed by: Matthew Bone
# Last Updated: 30/07/2021
# Updated by: Matthew Bone
#
# Contact Details:
# Bristol Composites Institute (BCI)
# Department of Aerospace Engineering - University of Bristol
# Queen's Building - University Walk
# Bristol, BS8 1TR
# U.K.
# Email - matthew.bone@bristol.ac.uk
#
# File Description:
# A unit test file designed for PyTest. This tests the neighbour searching tools
##############################################################################

from LammpsSearchFuncs import get_neighbours, get_additional_neighbours

atomList = ['1', '2', '3', '4', '5', '6', '7', '8']
bondList = [['1', '1', '1', '2'], ['2', '1', '2', '3'], ['3', '1', '2', '4'], ['4', '1', '2', '5'], ['5', '1', '1', '6'], ['6', '1', '6', '7'], ['7', '1', '7', '8']]
bondingAtoms = ['1', '2']

# Pseudochemistry for neighbours is:
# 8 - 7 - 6 - 1 - 2 - 3/4/5 # Think 8761 as carbon chain and 2345 as a methyl group

def inspect_values(result, expected):
    outputList = []
    for value in result:
        if value in expected:
            outputList.append(True)
        else: 
            outputList.append(False)

    return outputList

def test_get_neighbours():
    neighboursDict = get_neighbours(atomList, bondList, bondingAtoms)

    symmetryCheck = neighboursDict['3'] == neighboursDict['4'] and neighboursDict['4'] == neighboursDict['5']

    checkValues = [neighboursDict['1'], neighboursDict['2'], neighboursDict['3'], symmetryCheck, neighboursDict['7']] 
    expected = [['2', '6'], ['1', '3', '4', '5'], ['2'], True, ['6', '8']]

    assert checkValues == expected

def test_get_second_neighbours():
    # Inital requirements
    neighboursDict = get_neighbours(atomList, bondList, bondingAtoms)
    neighboursOne = neighboursDict[atomList[0]]
    neighboursTwo = neighboursDict[atomList[1]]
    neighboursSeven = neighboursDict[atomList[6]]

    secondNeighboursOne = get_additional_neighbours(neighboursDict, atomList[0], neighboursOne, bondingAtoms)
    secondNeighboursTwo = get_additional_neighbours(neighboursDict, atomList[1], neighboursTwo, bondingAtoms)
    secondNeighboursSevenAll = get_additional_neighbours(neighboursDict, atomList[6], neighboursSeven, bondingAtoms, unique=False)

    oneCheck = inspect_values(secondNeighboursOne, ['3', '4', '5', '7'])
    twoCheck = inspect_values(secondNeighboursTwo, ['6'])
    sevenCheckAll = inspect_values(secondNeighboursSevenAll, ['7', '1'])

    checkValues = [all(oneCheck), all(twoCheck), all(sevenCheckAll)] 
    expected = [True, True, True]

    assert checkValues == expected

def test_get_third_neighbours():
    # Inital requirements
    neighboursDict = get_neighbours(atomList, bondList, bondingAtoms)
    neighboursOne = neighboursDict[atomList[0]]
    secondNeighboursOne = get_additional_neighbours(neighboursDict, atomList[0], neighboursOne, bondingAtoms)
    
    thirdNeighboursOne = get_additional_neighbours(neighboursDict, atomList[0], searchNeighbours=secondNeighboursOne, bondingAtoms=bondingAtoms)
    thirdNeighboursOneAll = get_additional_neighbours(neighboursDict, atomList[0], searchNeighbours=secondNeighboursOne, bondingAtoms=bondingAtoms, unique=False)

    oneCheck = inspect_values(thirdNeighboursOne, ['8'])
    oneCheckAll = inspect_values(thirdNeighboursOneAll, ['8', '6', '2'])

    checkValues = [all(oneCheck), all(oneCheckAll)] 
    expected = [True, True]

    assert checkValues == expected