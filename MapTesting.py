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

import os
from PathSearch import map_from_path

class Reaction:
    def __init__(self, directory, preFileName, postFileName, elementByType):
        startDir = os.getcwd() # Save primary main directory path to return to it after mapping
        self.mappedIDList = map_from_path(directory, preFileName, postFileName, elementByType, True, mappedIDList=[])
        os.chdir(startDir)

    def test_report(self, correctPostAtomIDs, reactionName):
        print(f'Reaction: {reactionName}')
        # Print test report
        for mappedPair in self.mappedIDList:
            print(f'Atom {mappedPair[0]} is mapped to atom {mappedPair[1]}')

        
        totalAtoms = len(correctPostAtomIDs)
        correctAtoms = 0
        incorrectPreAtomsList = []
        for atom in self.mappedIDList:
            if atom[1] in correctPostAtomIDs[atom[0]]:
                correctAtoms += 1
            else:
                incorrectPreAtomsList.append(atom[0])

        mappedPostAtomsList = [val[1] for val in self.mappedIDList]
        repeatedPostIDs = [val for val in mappedPostAtomsList if mappedPostAtomsList.count(val) > 1]

        print(f'Total atoms: {totalAtoms}. Correct atoms: {correctAtoms}. Accuracy: {round(correctAtoms / totalAtoms * 100, 1)}%')
        print(f'Incorrectly assigned premolecule atomIDs: {incorrectPreAtomsList}, Count {len(incorrectPreAtomsList)}')
        print(f'Repeated Atoms: {repeatedPostIDs}, Count: {len(repeatedPostIDs)}\n\n')

# DGEBA-DETDA
dgebaDetda = Reaction(
    'Test_Cases/Map_Tests/DGEBA_DETDA/', 
    'pre-molecule.data', 
    'post-molecule.data', 
    ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O']
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
dgebaDetda.test_report(correctDgebaDetda, 'DGEBA-DETDA')

# Ethyl Ethanoate
ethylEthanoate = Reaction(
    'Test_Cases/Map_Tests/Ethyl_Ethanoate', 
    'pre-molecule.data', 
    'post-molecule.data', 
    ['H', 'H', 'C', 'C', 'O', 'O', 'O', 'O']
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
ethylEthanoate.test_report(correctEthylEthanoate, 'Ethyl Ethanoate')

# Methane to Ethane 
ethane = Reaction(
    'Test_Cases/Map_Tests/Methane_Ethane', 
    'pre-molecule.data', 
    'post-molecule.data', 
    ['H', 'C']
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
ethane.test_report(correctEthane, 'Methane to Ethane')

# LAMMPS Example - Nylon 6,6 taken from 'nylon,6-6_melt' example
nylon = Reaction(
    'Test_Cases/Map_Tests/Lammps_Nylon', 
    'rxn1_stp1_unreacted.data_template', 
    'rxn1_stp1_reacted.data_template', 
    ['C', 'N', 'H', 'H', 'C', 'O', 'H', 'O', 'N', 'H', 'O']
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
nylon.test_report(correctNylon, 'Nylon Melt Lammps Example')

# Phenol O-Alkylation
phenAlkyl = Reaction(
    'Test_Cases/Map_Tests/Phenol_Alkylation', 
    'pre-molecule.data', 
    'post-molecule.data', 
    ['H', 'H', 'C', 'C', 'O', 'O']
)
correctPhenAlkyl = {
    '1': ['1'],
    '2': ['2'],
    '3': ['3'],
    '4': ['4'],
    '5': ['5'],
    '6': ['6'],
    '7': ['7'],
    '8': ['8'],
    '9': ['9'],
    '10': ['10'],
    '11': ['11'],
    '12': ['23', '24'],
    '13': ['13'],
    '14': ['14'],
    '15': ['17'],
    '16': ['12', '15', '16'],
    '17': ['12', '15', '16'],
    '18': ['12', '15', '16'],
    '19': ['18', '19', '23', '24'],
    '20': ['18', '19', '23', '24'],
    '21': ['18', '19', '23', '24'],
    '22': ['20'],
    '23': ['21', '22'],
    '24': ['21', '22']
}
phenAlkyl.test_report(correctPhenAlkyl, 'Phenol O-Alkylation')

# Symmetric Diol
symmDiol = Reaction(
    'Test_Cases/Map_Tests/Symmetric_Diol', 
    'pre-molecule.data', 
    'post-molecule.data', 
    ['H', 'H', 'C', 'C', 'O', 'O']
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
symmDiol.test_report(correctSymmDiol, 'Symmetric Diol')
