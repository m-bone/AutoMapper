from PathSearch import map_from_path

class Reaction:
    def __init__(self, directory, preFileName, postFileName, elementByType):
        self.mappedIDList = map_from_path(directory, preFileName, postFileName, elementByType)

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
    '/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/FullModel/Reaction', 
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
    '/home/matt/Documents/Oct20-Dec20/Bonding_Test/Ethyl_Ethanoate/Reaction', 
    'pre-molecule.data', 
    'post-molecule.data', 
    ['H', 'H', 'C', 'C', 'O', 'O', 'O', 'O']
)
correctEthylEthanoate = {
    '11': ['2'],
    '6': ['7'],
    '10': ['1'],
    '12': ['3', '4', '5'],
    '13': ['3', '4', '5'],
    '14': ['3', '4', '5'],
    '17': ['6'],
    '2': ['8'],
    '7': ['10', '11'],
    '8': ['10', '11'],
    '1': ['9'],
    '3': ['12', '13', '14'],
    '4': ['12', '13', '14'],
    '5': ['12', '13', '14'],
    # Water molecule atoms
    '9': ['17', '16'],
    '16': ['17', '16'],
    '15': ['15']
}
ethylEthanoate.test_report(correctEthylEthanoate, 'Ethyl Ethanoate')

# 15pre maps to 6post instead of 17pre. This is correct as far as the code is concerned, but is the wrong map.
# Need to test what LAMMPS does when it receives a wrong assignment like this - will it error?


# Methane to Ethane 
ethane = Reaction(
    '/home/matt/Documents/Oct20-Dec20/Bonding_Test/Methane_Ethane/Reaction', 
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

# LAMMPS Example - Nylon 6,6
nylon = Reaction(
    '/home/matt/Documents/Oct20-Dec20/Bonding_Test/Nylon6-6', 
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
    '/home/matt/Documents/Oct20-Dec20/Bonding_Test/Phenol_Alkylation/Reaction', 
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
    '19': ['23', '24'],
    '20': ['18', '19'],
    '21': ['18', '19'],
    '22': ['20'],
    '23': ['21', '22'],
    '24': ['21', '22']
}
phenAlkyl.test_report(correctPhenAlkyl, 'Phenol O-Alkylation')

# Validation idea
# Search for ambiguous groups in the post molecule by comparing post atom to all post atoms
# I can find this easily and it can confirm if something should be an ambiguous group
# Could cause issues if the BPDM manages to split two things that should be pairs - this check may find things that the BPDM doesn't
# Tool would help explin why BPDM works in some cases but less in others
# I can also predict how many ambiguous groups in my pre and post molecule with this method
# This could check if I have as many as I expect and it may help identify atoms that have moved - can use ambiguous pairs as a useful tool
# Can I use ambiguous pairs that don't exist before but do after and visa versa to identify moved atoms