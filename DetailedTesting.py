from AtomMapping import atom_mapping

class Reaction:
    def __init__(self, directory, preFileName, postFileName, elementByType, preBondingAtoms, postBondingAtoms, postMajorMovedAtoms, postMinorMovedAtoms):
        self.mappedIDList = atom_mapping(directory, preFileName, postFileName, elementByType, preBondingAtoms, postBondingAtoms, postMajorMovedAtoms, postMinorMovedAtoms)

    def test_report(self, correctPostAtomIDs, reactionName):
        print(f'\n\nReaction: {reactionName}')
        # Print test report
        for mappedPair in self.mappedIDList:
            print(f'Atom {mappedPair[0]} is mapped to atom {mappedPair[1]}')

        
        totalAtoms = len(correctPostAtomIDs)
        correctAtoms = 0
        incorrectPreAtomsList = []
        for index, atom in enumerate(self.mappedIDList):
            if atom[1] in correctPostAtomIDs[index]:
                correctAtoms += 1
            else:
                incorrectPreAtomsList.append(atom[0])

        mappedPostAtomsList = [val[1] for val in self.mappedIDList]
        repeatedPostIDs = [val for val in mappedPostAtomsList if mappedPostAtomsList.count(val) > 1]

        print(f'Total atoms: {totalAtoms}. Correct atoms: {correctAtoms}. Accuracy: {round(correctAtoms / totalAtoms * 100, 1)}%')
        print(f'Incorrect premolecule atomIDs: {incorrectPreAtomsList}')
        print(f'Repeated Atoms: {repeatedPostIDs}, Count: {len(repeatedPostIDs)}')

# DGEBA-DETDA
dgebaDetda = Reaction('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction', 'new_start_molecule.data', 'new_post_rx1_molecule.data', ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O'],
['28', '62'], ['32', '15'], ['33'], ['16'])
correctDgebaDetda = [['38'], ['39'], ['35'], ['41', '42'], ['42', '41'], ['32'], ['16'], ['5', '36'], ['36', '5'], ['37'], ['6', '9'], ['4'], ['1', '3'], ['3', '1'], ['9', '6'], ['17', '23'], ['23', '17'], ['15'], ['33', '34'], ['34', '33']]
dgebaDetda.test_report(correctDgebaDetda, 'DGEBA-DETDA')

# Ethyl Ethanoate
ethylEthanoate = Reaction('/home/matt/Documents/Oct20-Dec20/Bonding_Test/Ethyl_Ethanoate/Reaction', 'pre-molecule.data', 'post-molecule.data', ['H', 'H', 'C', 'C', 'O', 'O', 'O', 'O'], ['6', '11'], ['7', '2'], ['17', '15', '16'], [])
correctEthylEthanoate = [['9'], ['8'], ['12', '13', '14'], ['13', '12', '14'], ['14', '12', '13'], ['7'], ['10', '11'], ['11', '10'], ['17', '16'], ['1'], ['2'], ['3', '4', '5'], ['4', '3', '5'], ['5', '3', '4'], ['15'], ['16', '17'], ['6']]
ethylEthanoate.test_report(correctEthylEthanoate, 'Ethyl Ethanoate')

# Nothing reasonable got given 13, too many 14 including some across the molecule boundary
# 15 given 2 should be impossible for multiple reasons - 15 is O, 2 is C and 2 is a bonding atom

# Validation idea
# Search for ambiguous groups in the post molecule by comparing post atom to all post atoms
# I can find this easily and it can confirm if something should be an ambiguous group
# Could cause issues if the BPDM manages to split two things that should be pairs - this check may find things that the BPDM doesn't
# Tool would help explin why BPDM works in some cases but less in others
# I can also predict how many ambiguous groups in my pre and post molecule with this method
# This could check if I have as many as I expect and it may help identify atoms that have moved - can use ambiguous pairs as a useful tool
# Can I use ambiguous pairs that don't exist before but do after and visa versa to identify moved atoms