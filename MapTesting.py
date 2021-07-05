##############################################################################
# Developed by: Matthew Bone
# Last Updated: 16/06/2021
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
# is given. These tests are not designed to be chemically practical; they're
# designed to tax different parts of the path search system to check functionality
##############################################################################

from PathSearch import map_from_path
from MapProcessor import map_processor, restore_dir

# Toggle Debug Reports
DEBUG = True

def test_report(mappedIDList, correctPostAtomIDs, reactionName):
    print(f'Reaction: {reactionName}')
    # Print test report
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
        ['28', '65'], None, ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O'], debug=DEBUG
    )

correctDgebaDetda = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '4': ['4'],
    '5': ['18'],
    '6': ['5'],
    '7': ['6'],
    '8': ['7'],
    '9': ['8'],
    '10': ['9'],
    '11': ['10'],
    '12': ['11'],
    '13': ['12', '15'],
    '14': ['13'],
    '15': ['14'],
    '16': ['15', '12'],
    '17': ['19', '20'],
    '18': ['16'],
    '19': ['20', '19'],
    '20': ['17']
}
test_report(ddMappedIDList, correctDgebaDetda, 'DGEBA-DETDA')

# Ethyl Ethanoate
with restore_dir():
    eeMappedIDList, _, _, _, _ = map_from_path(
        'Test_Cases/Map_Tests/Ethyl_Ethanoate/', 'pre-molecule.data', 'post-molecule.data', ['H', 'H', 'C', 'C', 'O', 'O', 'O', 'O'], debug=DEBUG
    )
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
test_report(eeMappedIDList, correctEthylEthanoate, 'Ethyl Ethanoate')

# Methane to Ethane
with restore_dir():
    meMappedIDList, _, _, _, _ = map_from_path(
        'Test_Cases/Map_Tests/Methane_Ethane/', 'pre-molecule.data', 'post-molecule.data', ['H', 'C'], debug=DEBUG
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
test_report(meMappedIDList, correctEthane, 'Methane to Ethane')

# LAMMPS Example - Nylon 6,6 taken from 'nylon,6-6_melt' example
with restore_dir():
    lnMappedIDList, _, _, _, _ = map_from_path(
        'Test_Cases/Map_Tests/Lammps_Nylon/', 'rxn1_stp1_unreacted.data_template', 'rxn1_stp1_reacted.data_template', ['C', 'N', 'H', 'H', 'C', 'O', 'H', 'O', 'N', 'H', 'O'], debug=DEBUG
    )

correctNylon = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '4': ['4', '5'],
    '5': ['4', '5'],
    '6': ['6', '7'],
    '7': ['6', '7'],
    '8': ['8'],
    '9': ['9'],
    '10': ['10'],
    '11': ['11'],
    '12': ['12'],
    '13': ['13', '14'],
    '14': ['13', '14'],
    '15': ['15'],
    '16': ['16'],
    '17': ['17', '18'],
    '18': ['17', '18'],
}
test_report(lnMappedIDList, correctNylon, 'Nylon Melt Lammps Example')

# Phenol O-Alkylation
with restore_dir():
    paMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Phenol_Alkylation/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['13', '14'], 
        ['13', '14'], ['12', '19', '23', '24'], ['H', 'H', 'C', 'C', 'O', 'O'], debug=DEBUG
    )

correctPhenAlkyl = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '4': ['4'],
    '5': ['5'],
    '6': ['6'],
    '7': ['7'],
    '8': ['19', '20'], # O-H hydrogen
    '9': ['9'],
    '10': ['10'],
    '11': ['13'],
    '12': ['8', '11', '12'],
    '13': ['8', '11', '12'],
    '14': ['8', '11', '12'],
    '15': ['20', '19', '14', '15'], # Bonding carbon H
    '16': ['20', '19', '14', '15'], # Bonding carbon H
    '17': ['20', '19', '14', '15'], # Bonding carbon H
    '18': ['16'],
    '19': ['18', '17'],
    '20': ['18', '17'],
}
test_report(paMappedIDList, correctPhenAlkyl, 'Phenol O-Alkylation')

# Symmetric Diol
with restore_dir():
    sdMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Symmetric_Diol/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '16'], 
        ['1', '16'], None, ['H', 'H', 'C', 'C', 'O', 'O'], debug=DEBUG
    )

correctSymmDiol = {
    '1': ['1'],
    '2': ['15'],
    '3': ['3'],
    '4': ['4', '14', '16', '2'],
    '5': ['5'],
    '6': ['6', '12'],
    '7': ['7', '11'],
    '8': ['8'],
    '9': ['9', '10'],
    '10': ['9', '10'],
    '11': ['7', '11'],
    '12': ['6', '12'],
    '13': ['13'],
    '14': ['4', '14', '16', '2'],
    '15': ['4', '14', '16', '2'],
    '16': ['4', '14', '16', '2'],
    '17': ['17']
}
test_report(sdMappedIDList, correctSymmDiol, 'Symmetric Diol')

# Generic PU
with restore_dir():
    gpMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Generic_PU/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '36'], 
        ['1', '36'], None, ['H', 'H', 'C', 'C', 'C', 'C', 'N', 'N', 'O', 'O', 'O', 'O'], debug=DEBUG
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
    '9': ['9'],
    '10': ['10'],
    '11': ['11'],
    '12': ['12'],
    '13': ['13', '17'],
    '14': ['14'],
    '15': ['19'],
    '16': ['16'],
    '17': ['17', '13'],
    '18': ['15', '18'],
    '19': ['15', '18']
}
test_report(gpMappedIDList, correctGenPU, 'Generic PU')

# Edge Atom Symmetry
# The key test for this is that 8 and 9, and 11 and 12 are not assigned by inference, but with edge atom symmetry
with restore_dir():
    eaMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Edge_Atom_Symmetry/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '32'], 
        ['1', '32'], None, ['H', 'H', 'C', 'C', 'O', 'O'], debug=DEBUG
    )

correctEdgSym = {
    '1': ['1'],
    '2': ['17'],
    '3': ['3'],
    '4': ['19'],
    '5': ['5', '6'],
    '6': ['6', '5'],
    '7': ['7'],
    '8': ['8'],
    '9': ['9'],
    '10': ['10'],
    '11': ['11'],
    '12': ['12'],
    '13': ['13'],
    '14': ['14'],
    '15': ['15'],
    '16': ['16'],
    '17': ['18', '4', '2'],
    '18': ['18', '4', '2'],
    '19': ['18', '4', '2']
}
test_report(eaMappedIDList, correctEdgSym, 'Edge Atom Symmetry')

# Queue Tester
# If this were to do the edge atoms too late, it would have to infer the symmetry atoms 5, 7 and 8
with restore_dir():
    qtMappedIDList = map_processor(
        'Test_Cases/Map_Tests/Queue_Tester/', 'cleanedpre_reaction.data', 'cleanedpost_reaction.data', 'pre-molecule.data', 'post-molecule.data', ['1', '33'], 
        ['1', '33'], None, ['H', 'H', 'C', 'C', 'O', 'O'], debug=DEBUG
    )

correctQTest = {
    '1': ['1'],
    '2': ['20'],
    '3': ['3'],
    '4': ['22'],
    '5': ['5'],
    '6': ['6', '13'],
    '7': ['7'],
    '8': ['8'],
    '9': ['9'],
    '10': ['10', '11'],
    '11': ['11', '10'],
    '12': ['12'],
    '13': ['13', '6'],
    '14': ['14'],
    '15': ['15', '16'],
    '16': ['16', '15'],
    '17': ['17'],
    '18': ['18'],
    '19': ['19'],
    '20': ['21', '4', '2'],
    '21': ['21', '4', '2'],
    '22': ['21', '4', '2'],

}
test_report(qtMappedIDList, correctQTest, 'Edge Atom Symmetry')