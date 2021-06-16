import os
from LammpsUnifiedCleaner import file_unifier
from LammpsTreatmentFuncs import clean_data
from LammpsSearchFuncs import get_data, find_sections

def test_unified_cleaner():
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Test_Cases/Cleaner/Methane_Ethane/') # Allows for relative pathing in pytest
    file_unifier(path, 'system.in.settings', ['pre-system.data', 'post-system.data'])

    # Load cleaned pre-system
    with open('cleanedpre-system.data', 'r') as f:
        data = f.readlines()

    data = clean_data(data)
    sectionIndex = find_sections(data)
    atoms = get_data('Atoms', data, sectionIndex)
    matchAtomsCount = int(data[1].split()[0]) == len(atoms)

    # Load cleaned settings
    with open('cleanedsystem.in.settings', 'r') as f:
        settings = f.readlines()

    # Number of coeffs, number of sections (+1 for end of file), last section name, number of bond coeffs, number of atoms in header and number in Atoms section
    checkValues = [len(settings), len(sectionIndex), data[sectionIndex[-2]], int(data[7][0]), matchAtomsCount] 
    expected = [8, 5, 'Angles', 3, True]

    assert checkValues == expected