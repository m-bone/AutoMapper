import os
from LammpsSearchFuncs import get_neighbours, get_additional_neighbours, get_data, find_partial_structure, find_sections
from LammpsTreatmentFuncs import clean_data

atomList = ['1', '2', '3', '4', '5', '6', '7', '8']
bondList = [['1', '1', '1', '2'], ['2', '1', '2', '3'], ['3', '1', '2', '4'], ['4', '1', '2', '5'], ['5', '1', '1', '6'], ['6', '1', '6', '7'], ['7', '1', '7', '8']]

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
    neighboursDict = get_neighbours(atomList, bondList)

    symmetryCheck = neighboursDict['3'] == neighboursDict['4'] and neighboursDict['4'] == neighboursDict['5']

    checkValues = [neighboursDict['1'], neighboursDict['2'], neighboursDict['3'], symmetryCheck, neighboursDict['7']] 
    expected = [['2', '6'], ['1', '3', '4', '5'], ['2'], True, ['6', '8']]

    assert checkValues == expected

def test_get_second_neighbours():
    # Inital requirements
    neighboursDict = get_neighbours(atomList, bondList)
    neighboursOne = neighboursDict[atomList[0]]
    neighboursTwo = neighboursDict[atomList[1]]
    neighboursSeven = neighboursDict[atomList[6]]

    secondNeighboursOne = get_additional_neighbours(neighboursDict, atomList[0], neighboursOne)
    secondNeighboursTwo = get_additional_neighbours(neighboursDict, atomList[1], neighboursTwo)
    secondNeighboursSevenAll = get_additional_neighbours(neighboursDict, atomList[6], neighboursSeven, unique=False)

    oneCheck = inspect_values(secondNeighboursOne, ['3', '4', '5', '7'])
    twoCheck = inspect_values(secondNeighboursTwo, ['6'])
    sevenCheckAll = inspect_values(secondNeighboursSevenAll, ['7', '1'])

    checkValues = [all(oneCheck), all(twoCheck), all(sevenCheckAll)] 
    expected = [True, True, True]

    assert checkValues == expected

def test_get_third_neighbours():
    # Inital requirements
    neighboursDict = get_neighbours(atomList, bondList)
    neighboursOne = neighboursDict[atomList[0]]
    secondNeighboursOne = get_additional_neighbours(neighboursDict, atomList[0], neighboursOne)
    
    thirdNeighboursOne = get_additional_neighbours(neighboursDict, atomList[0], searchNeighbours=secondNeighboursOne)
    thirdNeighboursOneAll = get_additional_neighbours(neighboursDict, atomList[0], searchNeighbours=secondNeighboursOne, unique=False)

    oneCheck = inspect_values(thirdNeighboursOne, ['8'])
    oneCheckAll = inspect_values(thirdNeighboursOneAll, ['8', '6', '2'])

    checkValues = [all(oneCheck), all(oneCheckAll)] 
    expected = [True, True]

    assert checkValues == expected

def test_get_partial_structure():
    os.chdir(path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Test_Cases/DGEBA_DETDA/')) # Allows for relative pathing in pytest

    # Load file into python as a list of lists
    with open('cleanedpre-reaction.data', 'r') as f:
        lines = f.readlines()
    
    # Prepare input
    tidiedLines = clean_data(lines)
    sectionIndexList = find_sections(tidiedLines)
    originalBonds = get_data('Bonds', tidiedLines, sectionIndexList)
    bondingAtoms = ['65', '28'] # 'C', 'N'

    validAtomSet, edgeAtomList, edgeAtomFingerprintDict = find_partial_structure(bondingAtoms, originalBonds, bondDistance=3)

    atomCheck = inspect_values(validAtomSet, ['1', '2', '3', '6', '8', '9', '13', '16', '28', '29', '37', '63', '64', '65', '66', '67', '68', '69', '70', '71'])
    edgeAtomCheck = inspect_values(edgeAtomList, ['1', '3', '13', '16', '37'])
    oneEdgeFingerprint = inspect_values(edgeAtomFingerprintDict['1'], ['4', '30'])
    thirteenEdgeFingerprint = inspect_values(edgeAtomFingerprintDict['13'], ['10', '14', '15'])

    checkValues = [len(validAtomSet), all(atomCheck), all(edgeAtomCheck), all(oneEdgeFingerprint), all(thirteenEdgeFingerprint)] 
    expected = [20, True, True, True, True]

    assert checkValues == expected

