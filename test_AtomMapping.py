from AtomMapping import atom_mapping

def validation_function(mappedIDList, correctPostAtomIDs):
    # Calculate accuracy
    totalAtoms = len(correctPostAtomIDs)
    correctAtoms = 0
    incorrectPreAtomsList = []
    for index, atom in enumerate(mappedIDList):
        if atom[1] in correctPostAtomIDs[index]:
            correctAtoms += 1
        else:
            incorrectPreAtomsList.append(atom[0])

    accuracy = round(correctAtoms / totalAtoms * 100, 1)

    # Calculate multiple assignment atoms
    mappedPostAtomsList = [val[1] for val in mappedIDList]
    repeatedPostIDs = [val for val in mappedPostAtomsList if mappedPostAtomsList.count(val) > 1]
    countRepeatedPostIDs = len(repeatedPostIDs)

    return accuracy, countRepeatedPostIDs

def test_dgeba_detda():
    mappedIDList = atom_mapping('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction', 'new_start_molecule.data', 'new_post_rx1_molecule.data', ['H', 'H', 'C', 'C', 'N', 'O', 'O', 'O'],
    ['28', '62'], ['32', '15'], ['33'], ['16'])
    correctPostAtomIDs = [['38'], ['39'], ['35'], ['41', '42'], ['42', '41'], ['32'], ['16'], ['5', '36'], ['36', '5'], ['37'], ['6', '9'], ['4'], ['1', '3'], ['3', '1'], ['9', '6'], ['17', '23'], ['23', '17'], ['15'], ['33', '34'], ['34', '33']]
    acc, repeatCount = validation_function(mappedIDList, correctPostAtomIDs)

    # Check accuracy and number of repeated IDs are as expect
    checkValues = [acc, repeatCount] 
    expected = [95, 2]

    assert checkValues == expected