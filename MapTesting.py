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
# Runs a series of different test molecules and print out resultant bonding pairs.
# These pairs are compared to the known correct pairs and an accuracy percentage
# is given. Many of these tests are not designed to be chemically practical; they're
# designed to tax different parts of the path search system to check functionality
##############################################################################

from MapProcessor import map_processor, restore_dir

# Toggle Debug and Test Reports
DEBUG = False
TEST_DEBUG = False

def test_report(mappedIDList, correctPostAtomIDs, reactionName, reactionForm):
    print(f'Reaction: {reactionName}')
    if reactionForm == 'Full':
        mappedIDList = mappedIDList[0]
    elif reactionForm == 'Partial':
        mappedIDList = mappedIDList[1]

    # Print test report
    if TEST_DEBUG:
        for mappedPair in mappedIDList:
            print(f'Atom {mappedPair[0]} is mapped to atom {mappedPair[1]}')

    
    totalAtoms = len(correctPostAtomIDs)
    correctAtoms = 0
    incorrectPreAtomsList = []
    for atom in mappedIDList:
        if atom[1] in correctPostAtomIDs[atom[0]]:
            correctAtoms += 1
        else:
            incorrectPreAtomsList.append(atom[0])

    mappedPostAtomsList = [val[1] for val in mappedIDList]
    repeatedPostIDs = [val for val in mappedPostAtomsList if mappedPostAtomsList.count(val) > 1]

    print(f'Total atoms: {totalAtoms}. Correct atoms: {correctAtoms}. Accuracy: {round(correctAtoms / totalAtoms * 100, 1)}%')
    print(f'Incorrectly assigned premolecule atomIDs: {incorrectPreAtomsList}, Count {len(incorrectPreAtomsList)}')
    print(f'Repeated Atoms: {repeatedPostIDs}, Count: {len(repeatedPostIDs)}\n\n')

# DGEBA-DETDA
with restore_dir():
    ddMappedIDList = map_processor(
        'Test_Cases/Map_Tests/DGEBA_DETDA/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['28', '65'], 
        ['28', '65'], None, ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O'], None, debug=DEBUG
    )

correctDgebaDetda = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '6': ['6'],
    '8': ['78', '29'],
    '9': ['9'],
    '13': ['13'],
    '16': ['16'],
    '28': ['28'],
    '29': ['78', '29'],
    '37': ['37'],
    '63': ['63'],
    '64': ['64', '67'],
    '65': ['65'],
    '66': ['66'],
    '67': ['64', '67'],
    '68': ['79', '80'],
    '69': ['69'],
    '70': ['79', '80'],
    '71': ['71'],
}
test_report(ddMappedIDList, correctDgebaDetda, 'DGEBA-DETDA', 'Partial')

# Ethyl Ethanoate
with restore_dir():
    eeMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Ethyl_Ethanoate/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['11', '6'],
        ['2', '7'], None, ['H', 'H', 'C', 'C', 'O', 'O', 'O', 'O', 'O', 'O'], None, debug=DEBUG
    ) # Del Atoms ['9', '15', '16', '16', '15', '17']
correctEthylEthanoate = {
    '1': ['9'],
    '2': ['8'],
    '3': ['12', '13', '14'],
    '4': ['12', '13', '14'],
    '5': ['12', '13', '14'],
    '6': ['7'],
    '7': ['10', '11'],
    '8': ['10', '11'],
    '10': ['1'],
    '11': ['2'],
    '12': ['3', '4', '5'],
    '13': ['3', '4', '5'],
    '14': ['3', '4', '5'],
    '17': ['6'],
    # Water molecule atoms
    '9': ['17', '16'],
    '16': ['17', '16'],
    '15': ['15']
}
test_report(eeMappedIDList, correctEthylEthanoate, 'Ethyl Ethanoate', 'Full')

# Methane to Ethane
with restore_dir():
    meMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Methane_Ethane/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '6'],
        ['1', '2'], ['5', '10', '9', '10'], ['H', 'C'], None, debug=DEBUG
    )
correctEthane = {
    '1': ['1'],
    '2': ['6', '7', '8'],
    '3': ['6', '7', '8'],
    '4': ['6', '7', '8'],
    '5': ['9', '10'],
    '6': ['2'],
    '7': ['3', '4', '5'],
    '8': ['3', '4', '5'],
    '9': ['3', '4', '5'],
    '10': ['9', '10'], 
}
test_report(meMappedIDList, correctEthane, 'Methane to Ethane', 'Full')

# Phenol O-Alkylation
with restore_dir():
    paMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Phenol_Alkylation/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['13', '14'], 
        ['13', '14'], ['12', '19', '23', '24'], ['H', 'H', 'C', 'C', 'O', 'O'], None, debug=DEBUG
    )

correctPhenAlkyl = {
    '1': ['1'],
    '2': ['2'],
    '4': ['4'],
    '5': ['5'],
    '6': ['6'],
    '9': ['9'],
    '10': ['10'],
    '12': ['23', '24'],
    '13': ['13'],
    '14': ['14'],
    '15': ['17'],
    '16': ['12', '15', '16'], 
    '17': ['12', '15', '16'],
    '18': ['12', '15', '16'],
    '19': ['23', '24'],
    '20': ['18', '19'],
    '21': ['18', '19'],
    '22': ['20'],
    '23': ['21', '22'],
    '24': ['21', '22'],
}
test_report(paMappedIDList, correctPhenAlkyl, 'Phenol O-Alkylation', 'Partial')

# Symmetric Diol
with restore_dir():
    sdMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Symmetric_Diol/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '16'], 
        ['1', '16'], None, ['H', 'H', 'C', 'C', 'O', 'O'], None, debug=DEBUG
    )

correctSymmDiol = {
    '1': ['1'],
    '2': ['18'],
    '5': ['5'],
    '7': ['7', '2', '17'],
    '8': ['8'],
    '9': ['9', '15'],
    '10': ['10'],
    '11': ['11'],
    '12': ['12', '13'],
    '13': ['12', '13'],
    '14': ['14'],
    '15': ['9', '15'],
    '16': ['16'],
    '17': ['7', '2', '17'],
    '18': ['7', '2', '17'],
    '19': ['19'],
    '20': ['20']
}
test_report(sdMappedIDList, correctSymmDiol, 'Symmetric Diol', 'Partial')

# Generic PU
with restore_dir():
    gpMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Generic_PU/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '36'], 
        ['1', '36'], None, ['H', 'H', 'C', 'C', 'C', 'C', 'N', 'N', 'O', 'O', 'O', 'O'], None, debug=DEBUG
    )

correctGenPU = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '4': ['4', '8'],
    '5': ['5'],
    '6': ['6'],
    '7': ['7'],
    '8': ['8', '4'],
    '11': ['11'],
    '12': ['12'],
    '30': ['30'],
    '31': ['31'],
    '32': ['32', '37'],
    '33': ['33'],
    '35': ['39'],
    '36': ['36'],
    '37': ['32', '37'],
    '38': ['35', '38'],
    '39': ['35', '38'],
}
test_report(gpMappedIDList, correctGenPU, 'Generic PU', 'Partial')

# Edge Atom Symmetry
# The key test for this is that 8 and 9, and 11 and 12 are not assigned by inference, but with edge atom symmetry
with restore_dir():
    eaMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Edge_Atom_Symmetry/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '32'], 
        ['1', '32'], None, ['H', 'H', 'C', 'C', 'O', 'O'], None, debug=DEBUG
    )

correctEdgSym = {
    '1': ['1'],
    '2': ['35'],
    '5': ['5'],
    '7': ['37'],
    '8': ['8', '11'],
    '11': ['8', '11'],
    '13': ['13'],
    '14': ['14'],
    '17': ['17'],
    '19': ['19'],
    '23': ['23'],
    '26': ['26'],
    '28': ['28'],
    '32': ['32'],
    '33': ['33'],
    '34': ['34'],
    '35': ['2', '7', '36'],
    '36': ['2', '7', '36'],
    '37': ['2', '7', '36'],
}
test_report(eaMappedIDList, correctEdgSym, 'Edge Atom Symmetry', 'Partial')

# Queue Tester
# If this were to do the edge atoms too late, it would have to infer the symmetry atoms 5, 7 and 8
with restore_dir():
    qtMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Queue_Tester/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '33'], 
        ['1', '33'], None, ['H', 'H', 'C', 'C', 'O', 'O'], None, debug=DEBUG
    )

correctQTest = {
    '1': ['1'],
    '2': ['36'],
    '5': ['5'],
    '7': ['38'],
    '8': ['8'],
    '9': ['9', '21'],
    '11': ['11'],
    '14': ['14'],
    '17': ['17'],
    '18': ['18', '19'],
    '19': ['18', '19'],
    '20': ['20'],
    '21': ['9', '21'],
    '27': ['27'],
    '28': ['28', '29'],
    '29': ['28', '29'],
    '33': ['33'],
    '34': ['34'],
    '35': ['35'],
    '36': ['2', '7', '37'],
    '37': ['2', '7', '37'],
    '38': ['2', '7', '37'],
}
test_report(qtMappedIDList, correctQTest, 'Queue Tester', 'Partial')

# Third Neighbour Symmetry
# Should determine atoms 8 and 11 by third neighbours, not inference
with restore_dir():
    tnMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Third_Neighbour_Symmetry/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', 
        ['1', '12'], ['1', '12'], None, ['H', 'H', 'C', 'C', 'O', 'O'], None, debug=DEBUG
    )

correctTNTest = {
    '1': ['1'],
    '2': ['33'],
    '3': ['3', '4'],
    '4': ['3', '4'],
    '5': ['5'],
    '7': ['35'],
    '8': ['8'],
    '11': ['11'],
    '12': ['12'],
    '13': ['13'],
    '14': ['14'],
    '15': ['2', '7', '34'],
    '16': ['2', '7', '34'],
    '17': ['2', '7', '34'],
    '18': ['18'],
    '23': ['23'],
    '26': ['26'],
    '28': ['28', '29'],
    '29': ['28', '29'],
}
test_report(tnMappedIDList, correctTNTest, 'Third Neighbour Symmetry', 'Partial')

# Caprolactam
with restore_dir():
    caMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Caprolactam/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', 
        ['3', '20'], ['3', '20'], None, ['H', 'H', 'C', 'C', 'N', 'N', 'N', 'N', 'O'], None, debug=DEBUG
    )

correctCATest = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '4': ['4', '5'],
    '5': ['4', '5'],
    '6': ['6', '7'],
    '7': ['6', '7'],
    '8': ['8'],
    '9': ['9'],
    '10': ['10', '11'],
    '11': ['10', '11'],
    '12': ['12'],
    '15': ['15'],
    '18': ['18'],
    '19': ['19'],
    '20': ['20'],
    '21': ['21'],
    '22': ['22'],
    '23': ['23', '24'],
    '24': ['23', '24'],
    '25': ['25', '26'],
    '26': ['25', '26'],
    '27': ['37'],
    '28': ['28'],
    '29': ['29', '30'],
    '30': ['29', '30'],
    '31': ['31'],
    '32': ['32', '33'],
    '33': ['32', '33'],
    '34': ['34'],
    '35': ['35', '36'],
    '36': ['35', '36'],
    '37': ['27'],
}
test_report(caMappedIDList, correctCATest, 'Caprolactam', 'Partial')

# Phenolic Resin
# This tests partial molecules with byproducts that aren't deleted
with restore_dir():
    prMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Phenolic_Resin/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', 
        ['4', '19'], ['4', '19'], None, ['H', 'H', 'C', 'C', 'O', 'O'], None, debug=DEBUG
    )

correctPRTest = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '4': ['4'],
    '6': ['6'],
    '7': ['7'],
    '8': ['8'],
    '9': ['9'],
    '10': ['10'],
    '13': ['29', '30'],
    '16': ['13', '28'],
    '17': ['13', '28'],
    '18': ['18'],
    '19': ['19'],
    '20': ['20'],
    '21': ['21'],
    '22': ['22'],
    '23': ['23'],
    '24': ['24'],
    '25': ['29', '30'],
    '26': ['26'],
    '27': ['27'],
    '29': ['17'],
    '30': ['25'],
}
test_report(prMappedIDList, correctPRTest, 'Phenolic Resin', 'Partial')

# Create Atoms Ethylene Glycol
# Tests how mapping handles atoms created in the post-bond structure.

with restore_dir():
    crMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Create_Atoms/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', 
        ['1', '2'], ['1', '2'], ['10', '25'], ['H', 'H', 'C', 'O', 'O'], createAtoms=['22', '10', '24', '20', '21', '18', '19', '23'],
        debug=DEBUG
    )

correctCRTest = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3', '8'],
    '4': ['4'],
    '5': ['5', '9'],
    '6': ['6', '15'],
    '7': ['7'],
    '8': ['8', '3'],
    '9': ['9', '5'],
    '10': ['25'],
    '11': ['11'],
    '14': ['14'],
    '15': ['15', '6'],
}

test_report(crMappedIDList, correctCRTest, 'Create Atoms', 'Partial')