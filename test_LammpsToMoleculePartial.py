from LammpsToMoleculePartial import lammps_to_molecule_partial
from LammpsTreatmentFuncs import clean_data
from LammpsSearchFuncs import get_data, find_sections, get_top_comments

def test_lammps_to_molecule_partial():
    lammps_to_molecule_partial('/home/matt/Documents/Bond_React_Python/Test_Cases', 'cleanedpost_rx1.data', 'post_rx1_', ['32', '15'])

    with open ('post_rx1_molecule.data', 'r') as f:
        mol = f.readlines()
    
    edgeAtoms = get_top_comments(mol)[1][1:]
    
    mol = clean_data(mol)
    sectionIndex = find_sections(mol)
    types = get_data('Types', mol, sectionIndex)
    charges = get_data('Charges', mol, sectionIndex)
    coords = get_data('Coords', mol, sectionIndex)
    atomSectionLengths = len(types) == len(charges) and len(types) == len(coords)

    # Number of section keywords (+1 for end of file), Last section header is Impropers, type of atom 6 is 2, atom based sections should be equal length, 12 angles in molecule
    # Number of edge atoms in the file
    checkValues = [len(sectionIndex), mol[sectionIndex[-2]], int(types[5][1]), atomSectionLengths, int(mol[2].split()[0]), len(edgeAtoms)] 
    expected = [8, 'Impropers', 4, True, 31, 5]

    assert checkValues == expected